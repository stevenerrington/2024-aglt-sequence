% % Subject 522:
% % ddir = fullfile('~/hbrl/HBRL_Analyses','Anal17_NLX-TDT_505-','Subject.522A.DBS.STN_TDT','Greenlee_Rohl');
%  ddir = fullfile('..','..','Anal17_NLX-TDT_505-','Subject.522A.DBS.STN_TDT','Greenlee_Rohl');
% % % File name is 522-002_SPECIALevents.mat. We will look at spike from Inpt_RZ2_chn001-Inpt_RZ2_chn003. The data is already high passed filtered.  Attached are the event time stamps for Speech production. 
% datafiles = fullfile(ddir,'522-002_SPECIALevents.mat');
% 
% % Subject 535:
% ddir = fullfile('..','..','Anal17_NLX-TDT_505-','Subject.535A.DBS.STN_TDT','Greenlee_Johari_OR_Timing');
% % file name  for STN is 535-001_SPECIALevents.mat , We will look at spikes from Inpt_RZ2_chn001-Inpt_RZ2_chn003. Input 3 is in STN
% datafiles = [datafiles,fullfile(ddir,{'535-001_SPECIALevents.mat','535-002_SPECIALevents.mat'})];
% % and for Nigra is 535-002_SPECIALevents.mat  we will look at Inpt_RZ2_chn001-Inpt_RZ2_chn003, 

% ddir = 'rawDemo'
ddir = 'osort-v4-rel/sim1'
dd = dir(fullfile(ddir,'simulation1.mat'));
datafiles = fullfile(ddir,{dd.name});
outputdir = 'hosd_output';

%%% load the data

%%% PARAMETERS FOR HOSD
windur = .008; % Window duration (sec)
lowpass = 4000 ; % lowpass cutoff ( Hz )
ncomp = 3; % Number of components to extract initially
window = 'hann'; %Window applied after segmentation of the input
tolerance = .001; % Tolerance for the purpose peak detection (width of smoothing kernel in sec. applied to peak detection output).
% skewness_threshold = .05; %Discard components whose filtered output has skewness below this threshold. 
skewness_z_threshold = 6; %Discard components for which the t-tstatistic computed on sample skewnesses falls below this value 
maxmem=16; %Max memory to use in GB. The spacing of window overlap samples is adjusted to stay within this value.

%%% INITIALIZE THE HOS OBJECT
% hos(ncomp) = hosobject(3); 
hos = hosobject(3); 
%%
for dfi = 1:length(datafiles)
    datafile = datafiles{dfi};
    
    for chan = 4
    %     chan = 1;
        hos = hos(1);
        hos.reset;
        [~,fn,ext] = fileparts(datafile);
        iterzthresh = 50;
        switch lower(ext) 
            case '.mat'
                if contains(fn,'simulation')
                    
                    ld = load(datafile);
                    vname = fn;
                    dat.dat = ld.spiketrains{chan};
                    dat.fs = 25000;
                    truspiket = [ld.spiketimes{:}];
                    truspikeval = arrayfun(@(x,k)ones(size(x{1}))*k,ld.spiketimes,1:3,'uniformoutput',false);
                else
                    vname = sprintf('Inpt_RZ2_chn%03i',chan);
                    dat = load(datafile,vname);
                    dat = dat.(vname);
                    iterzthresh = inf;
                end
            case '.ncs'
                dat = readncs(datafile);
                vname = fn;
                iterzthresh = 10;
        end
        hos.initialize(round(windur*dat.fs(1)),dat.fs(1),lowpass,[],[],'window',window,'glowpass',lowpass);

        %%% DENOISING: LINE NOISE AND TRANSIENTS
    %     figure
    %     [xdn,~,~,spk] = dbtDenoise(double(dat.dat),dat.fs(1),.25,'make plot',true);
        xdn = double(dat.dat);
        %%% Filter in MUA range
        gaps = convn(xdn~=0,ones(1000,1),'same')==0;
        xdn(gaps) =[];
        mua = hpfilt(zscore(xdn),[dat.fs(1) 300 ]);
    %     z = zscore(mua) + 0./(spk.filter>.5);
        z = iterz(zscore(mua),iterzthresh);
        
    %     z = z( t<=max(evt)+1);
        %%% APPLY HOSD TO DATA
        max_input_length = floor(maxmem*1e9/16/length(hos(1).B)/2*hos(1).buffersize/10);
        zrem = z;
        if length(z)>max_input_length
          warning('To stay within memory limits, sample window spacing will be increased by a factor of %0.2f',max(.5,length(zrem)./max_input_length)) 
        end
        go = true;
        keepc = 1;
        skewz=[];
        while go 

            fprintf('\nGetting component %i',keepc)

            hos(keepc) = hosobject(hos(1));
            if length(z)<=max_input_length
                [A,B] = hos(keepc).get_block(zrem); % A and B are segmented data after and before alignment, respectively.
            else
    %             for k = 1:length(hos)
                    %%% Deal with very long recordings by limiting the number of sample windows
                    hos(keepc).poverlap= max(.5,length(zrem)./max_input_length); %Accomodate the input size to memory limit
                                                                          %by adjusting the window spacing.
    %             end
                 [A,B] = hos(keepc).get_block(zrem); % A and B are segmented data after and before alignment, respectively.

        %         [A,B] = hos.get_block(z(1:max_input_length)); % Deal with very
        %                                                        long recordings by feeding in chunks of data sequentially
        %         for k = 1:length(hos)
        %             hos(k).burnin = 1e10;
        %         end
        %         for ti = max_input_length+1:max_input_length:length(z)
        %              hos.get_input(z(ti:min(end,ti+max_input_length-1)));
        %         end
            end
            %%% Apply the feature identification filter
            xfilt = hos(keepc).xfilt(zrem);
            xfilt(isnan(xfilt))=0;

            %%% function with which to compute skewness
            skewfun = @(x)nanmean(x.^3)./nanmean(x.^2).^(3/2);

        %     skw = skewfun(xfilt); % Compute skewness on filtered data    
        %     keepc = min(find(skw>=skewness_threshold,1,'last')+1,ncomp);
            %%% Use a statistical threshold on sample skewness to decide which
            %%% components to retain
            Tsamp = chopper(hos(1).segment.Trange,hos(1).segment.wint,hos(1).segment.fs);
            xskew =[];
            for k = 1:size(xfilt,2)
                xskew(:,k) = skewness(xfilt(Tsamp+(k-1)*size(xfilt,1)));
            end
            skewz(keepc) = nanmean(xskew)./nanstd(xskew)*sqrt(size(Tsamp,2));%Tstat on sample skewnes.
            zrem = zrem-hos(keepc).xrec(zrem);
            go = skewz(end)>skewness_z_threshold || keepc<3;
            keepc = keepc+go;
            close
        end
        keepc = find(skewz>skewness_z_threshold,1,'last');
        if isempty(keepc)
            keepc=1;
        end
        keepc = max(keepc,3);
    %     %%% Filtered and thresholded data
        xthr = hos(1:keepc).xthresh(z);

        %%%%% Get the peaks in the thresholded data %%%%%
        %%% Apply tolerance smoothing
        smxthr = convn(full(xthr),hann(round(2*tolerance*dat(1).fs(1))),'same');
        %%% Here combining over components and smoothing
        [~,pks] = getpeak2(sum(smxthr,2));
    %     [~,pks] = getpeak2(sum(xthr,2));
        %%% Here separately for each component
    %     [~,pksep] = getpeak2(smxthr); % With smoothing
        [~,pksep] = getpeak2(xthr); %No smoothing
        pks = pksep;
      
        pkts = [];
        compno = [];
        for k = 1:size(pks,2)
            pkts = [pkts,find(pks(:,k)'==1)];
            compno = [compno;k*ones(sum(pks(:,k)'==1),1)];
        end
        
        pksm = convn(full(pksep==1),ones(round(tolerance*dat.fs(1)),1),'same')>0;

        [pkts,unqi] = unique(pkts,'stable');
        compno = compno(unqi);
    %%%
        %compno = compno(diff([0,pkts])/dat.fs(1)>tolerance);
        
        [kni,knd] = knnsearch(pkt,pkt,'k',10);
        knD = sparse(length(pkt),length(pkt));
        knD(kni(:,2:end)+length(pkt)*(0:length(pkt)-1)') = knd(:,2:end)<tolerance*dat.fs(1);
        getkn = find(sum(knD)>0);
        G = graph(knD(getkn,getkn));
        L = laplacian(G);
        
        [V,D] = svd(full(L));
        
        
        %pkts = pkts(diff([0,pkts])/dat.fs(1)>tolerance);
        
    %     compno = pksm(pkts,:)==1;
    %     [~,compno] = max(zscore(smxthr(pkts,:)),[],2); 
    %     compno = compno.*(sum(smxthr(pkts,:),2)~=0);
        pkts = pkts(~all(compno==0,2));
    %     compno = compno(~all(compno==0,2),:);

        %%% Create segmentation matrix around spike times
        [T,tt] = chopper(windur*[-.5 .5]-[0 1/dat.fs(1)],pkts/dat.fs(1),dat.fs(1));
        pkts = pkts(:,~any(T<1 | T>length(z)));
        compno = compno(~any(T<1 | T>length(z)),:);
        T = T(:,~any(T<1 | T>length(z)));

        %%% Try clustering in the space of skewness  within the respective
        %%% peri-spike windows;
    %     X=[];
    %     for k = 1:ncomp
    %         xf = xfilt(:,k); 
    %         X(:,:,k) = xf(T);
    %     end

        if length(pkts)>sqrt(maxmem*1e9/8) %Limit number of spikes to a tractable value for clustering analyes
                                           %Component numbers still apply to
                                           %all.
            clusterpeaks =randperm(length(pkts)-1);
            clusterpeaks = clusterpeaks(1:floor(sqrt(maxmem*1e9/8)));
            clusterpeaks(end+1)=length(pkts);
        else
            clusterpeaks = 1:length(pkts);
        end


        %%% To reduce opportunity for bias related to different filter delays, realign waveforms to the first component according to the peak in the output of the 1st detection filter;
        [~,dt] = max(hos(1).xfilt([hos.feature]));
        dt = hos(1).sampt(dt(compno));
        Trl = T+dt;
        Trl(Trl<1)=1;Trl(Trl>length(z)) = length(z);

       xrecsep = hos(1:keepc).xrec(z);
       xrec = sum(xrecsep,2);
       xrfilt = hos.xfilt(xrec);
       xfilt = hos.xfilt(z);
        X = [];
        for k = 1:keepc
    %         X(:,:,k) = hos(k).xfilt(z(Trl).*hamming(size(T,1)));
             xrf = xfilt(:,k);
%              xrf = hos(k).xfilt(xrec);
    %         xrf = xthr(:,k);
             X(:,:,k) = xrf(T).*hamming(size(T,1));
    %         X(:,:,k) = hos(k).xfilt(xrec(T).*hamming(size(T,1)));
    %         xrfilt = hos(k).xfilt(xrecsep(:,k));
    %         X(:,:,k) = xrfilt(Trl).*hamming(size(T,1));
        end
        X(isnan(X))=0;
    %     Xsk = squeeze(skewfun(X));
    %     Xsk = squeeze(skewfun(X));

    %     Xsk = squeeze(max(X));
    %     Xsk = squeeze(skewness(X));

%         Xsk = zscore(squeeze(nthroot(sum(X(abs(tt)<.0005,:,:).^3),3)));

%         Xsk = zscore(squeeze(sum(X(abs(tt)<.0005,:,:).^3)));
        Xsk = zscore(squeeze(skewfun(X(abs(tt)<.001,:,:))));
%         Xsk = zscore(squeeze(sum(X(abs(tt)<.0001,:,:))));
    %     Xsk(:,end+1) = zscore(compno);
        if keepc==1
            Xsk=Xsk';
        end

    %     Xsk = Xsk(clusterpeaks,:,:);
        S = Xsk'*Xsk;
        [u3,l3] = svd(S);
        v3 = zscore(Xsk*u3*diag(1./sqrt(diag(l3))));

    %     XX = reshape(permute(zscore(X(abs(tt)<.0005,:,:)),[1 3 2]),sum(abs(tt)<.0005)*size(X,3),size(X,2));
    %     S2 = XX*XX';
    %     [u2,l2] = svd(S2);
    %     v2 = zscore(XX'*u2(:,1:ncomp)*diag(diag(l2(1:ncomp,1:ncomp)).^-.5));
    %     tic, cl3 = spectralcluster(zscore(v3),keepc);toc
    %     [kni,knd] = knnsearch(v3,v3,'K',50);

        %%%
    %     cluster_type = 'dbscan'
%         cluster_type = 'dbcan'
%         cluster_type = 'component';
        cluster_type = 'pca';
        switch cluster_type
            case 'pca'
                Zrl = z(Trl);
                hos2 = hosobject(hos);
                A = hos2.get_block(Zrl.*hann(size(T,1)));
                for kk = 1:size(A,3)
                    AF(:,:,kk) = zscore(fftshift(hos(kk).xfilt(A(:,:,kk)),1));
                end
                Ars = reshape(permute(AF(abs(tt)<.001,:,1:keepc),[1 3 2]),sum(abs(tt)<.001)*keepc,size(A,2));
                AA = (Ars*Ars');
                [u,l] = svd(AA);
                v = Ars'*u(:,1:keepc);
                kmfun = @(x,k)kmeans(x,k,'maxiter',500,'replicates',10);
                tic,eva = evalclusters(v3,kmfun,'CalinskiHarabasz','klist',2:15);toc
                cl3 = eva.OptimalY;
                for kk = 1:max(cl3)
                   hos3(kk,1) = hosobject(hos(1));
                   hos3(kk).get_block(z(T(:,cl3==kk))); 
                   Trl(:,cl3==kk)=T(:,cl3==kk)+hos3(kk).delay;
                end
                
            case 'kmeans'
                kmfun = @(x,k)kmeans(x,k,'maxiter',500,'replicates',10);
                tic,eva = evalclusters(v3,kmfun,'CalinskiHarabasz','klist',2:15);toc
                cl3 = eva.OptimalY;
                 cluster_var_threshold=Inf;
             case 'dbscan'
                scales = 2.^(1:.5:6);
                clust = zeros(size(v3,1),length(scales)+1);
                clust(:,1) = 1;
                hs = [];
                for k = 1:length(scales)
            %         dst = mean(knd(:,end))/scales(k);
                    dst = 1./scales(k);

                    if all(clust(:,k)==0)
                        break
                    end
                    [unq,~,unqi] = unique(clust,'rows');
                    for kk = 1:size(unq,1)
                        if unq(kk,k)==0
                            continue
                        end
                        vsub = v3(unqi==kk,:);
                        tic, cl3 = dbscan(vsub,dst,50);toc
                        cl3(cl3==-1)=0;
                        if ~any(cl3>0)
                            continue
                        end
                        [u2,~,u2i] = unique(cl3);
                        if max(u2)==1
                             clust(unqi==kk,k+1) = cl3;
            %                  clust(unqi==kk,k+1) = 1;
                             continue
                        end
                        mm=[];
                        S=[];
                        geti = find(u2>0)';
            %             geti = 1:length(u2);
                        for mk = 1:length(geti)
                            mm(mk,:) = mean(vsub(u2i==geti(mk),:));
                            S(:,:,mk) = cov(vsub(u2i==geti(mk),:))*sum(u2i==geti(mk));
                        end
                        mS = sum(S,3)./sum(cl3>0);
                        Sm = cov(mm);
                        h = svd(Sm*mS^-1);
                        hs(end+1)=h(1);
                        if h(1)<2
                            continue;
                        end
                        clust(unqi==kk,k+1) = cl3;

                    end
                    k
                end

                [unq,~,cl3] = unique(clust,'rows','stable');
                cluster_var_threshold=2;
            case 'umap'
                if isempty(which('umap'))
                    addpath UMAP/umap
                end
                [reduction, umap, clusterIdentifiers, extras]=run_umap(zscore(v3),'n_components',3,'min_dist',.1,'n_neighbors',50);
%                 [reduction, umap, clusterIdentifiers, extras]=run_umap(zscore(v3),'n_components',3,'min_dist',.3,'n_neighbors',min(199,round(sqrt(size(v3,1)))));
                cl3 = clusterIdentifiers;
                cluster_var_threshold=Inf;
            case 'component'
                cl3=compno ;
                cluster_var_threshold = Inf;
        end
        %%% Eliminate clusters whose centroids are separated by less than 1 stddev.
        run = true;
        while run

             clear SS MM S mhd
             [unq,~,cl3] = unique(cl3);

            for k = 1+(unq(1)==-1):max(cl3)
                S(:,:,k)= cov(v3(cl3==k,:));
                SS(k) = trace(cov(v3(cl3==k,:))); 
                MM(k,:) = mean(v3(cl3==k,:),1); 
            end

            for k = 1:size(MM,1)
                for kk = 1:size(MM,1)
                    mS = mean(S(:,:,[k kk]),3);
                    dm = diff(MM([k kk],:));
                    mhd(k,kk)  = dm*mS^-1*dm';
                end
            end        
    %         mhd=mhd(SS<1,SS<1);
    %          cl3(SS(cl3)>1)=-1; %Discard clusters with variance greater than 1
            [ii,jj] = ndgrid(1:size(mhd,1),1:size(mhd,1));
            ey = eye(size(mhd));
            [mn,mni]  = min(mhd(:) + 0./(1-ey(:)) + 0./(SS(ii(:))<cluster_var_threshold&SS(jj(:))<cluster_var_threshold)');
            [mi,mj] = ind2sub(size(mhd),mni);


            if mn<1 && sum(cl3==mi)>0
                cl3(cl3==mi) = mj;
            else
                run=false;
             end
    %         cl3(SS(cl3)>2)=-1; %Discard clusters with large variance 
    %         cl3(unq(1)==-1 & cl3==1)=-1;
           cl3(SS(cl3(:))>cluster_var_threshold | unq(cl3)==-1)=-1; %Discard clusters with large variance 
           cl3(cl3>0) = cl3(cl3>0)-(unq(1)==-1);
    %         mn
        end

        MHD = [];
        for k = 1:max(cl3)
            MHD(cl3==k) = mahal(v3(cl3==k,:),v3(cl3==k,:));
        end

    %      cl3(SS(cl3)>2)=-1; %Discard clusters with large variance 
    %      [unq,~,cl3] = unique(cl3);
    %     cl3(unq(cl3)==-1) = -1;
    %     cl3(cl3>0) = cl3(cl3>0)-(unq(1)==-1);
    %     klist=2:20;%the number of clusters you want to try
    %     
    %      myfunc = @(X,K)kmeans(X, K, 'emptyaction','singleton','replicate',10,'MaxIter',500);
    %     eva = evalclusters(zscore(v3),myfunc,'CalinskiHarabasz','klist',klist)
    %     cl3=eva.OptimalY;
    % 
    %     unq = unique(cl3);

        clear M3 M3x
         for k = 1:max(cl3)
    %          M3(:,k) = mean(Q(:,cl3==k),2);
    %          M3x(:,k) = mean(z(Trl(:,clusterpeaks(cl3==k))),2);
             M3x(:,k) = mean(z(Trl(:,cl3==k)),2);
         end

         clear R IMP
         for k = 1:max(cl3)

             q = zeros(round(size(z,1)*1000/dat.fs(1)),1);
             q(round(pkts(cl3==k)/dat.fs(1)*1e3))=1;
             IMP(:,k) = q;
             R(:,k) = convn(q,ones(30*1e3,1)/30,'same');

         end
          v3(:,end+1:3) = 0;  
        %%
        [srt,srti] = sort(cl3);

        fig1 = figure;
        subplot(1,3,1),
        cm = hsv(double(max(cl3)));
    %     cm = cm(randperm(size(cm,1)),:);
        plot3(v3(cl3<0,1),v3(cl3<0,2),v3(cl3<0,3),'.k','markersize',1)    
        xlabel PC1
        ylabel PC2
        zlabel PC3
        hold on
        scatter3(v3(cl3>0,1),v3(cl3>0,2),v3(cl3>0,3),1,cm(cl3(cl3>0),:))    
    %     colormap(cm)
        caxis(.5+[0 max(cl3)])
        axis vis3d
        grid on
        subplot(1,4,4)
    %     plot(tt,M3x+5*(1:size(M3,2)),'k')
        plot(tt([1 end])*1e3, [1 1]'*(1:max(srt(srt>0))),'k:')
        hold on 
        plh = plot(tt*1e3,M3x./(8*std(M3x(:)))+(1:size(M3x,2)),'r');
        xlim([-2 2])
        xlabel('ms')
        set(plh,'linewidth',2)
        title('Mean waveforms')
        subplot(2,3,2)
        imagesc(tt*1e3,[],z(Trl(:,srti(srt>0)))')
        hold on, plh(:,2)=plot(0./(srt(srt>0)==(1:max(srt(srt>0))))+2,1:sum(srt>0));
        title('Spike waveforms')
        xlim([-2 2])
        xlabel('ms')
        axis xy
        subplot(2,3,5)
        imagesc(tt*1e3,[],convn(z(Trl(:,srti(srt>0))),ones(1,100)/100,'same')')
        title('Averaged Waveforms (100 point moving avg)')
        xlabel('ms')
        xlim([-2 2])
        axis xy
        hold on, plh(:,3)=plot(0./(srt(srt>0)==(1:max(srt(srt>0))))+2,1:sum(srt>0));
        for k = 1:size(plh,1),set(plh(k,:),'Color',cm(k,:));end
    %     plh=plot(0./(srt==(1:max(srt)))+1);
        set(plh(:,2:3),'linewidth',5)
        set(fig1,'units','normalized','position',[0 0 1 1]);
        cl=[];
           %% 
    %     %%% Alternative: cluster according to waveforms after common alignment
    %     %%% This applies HOSD again on the waveforms and the waveforms aligned 
    %     %%% on the respective component
    % %     hos2 = hosobject(hos(1:keepc));
    % %     hos2(end+1) = hosobject(hos(1));
    % %     [A,B] = hos2.get_block(z(T));
    % 
    % %     Q = z(Trl).*hamming(size(T,1));
    %     Q = xrec(Trl(:,clusterpeaks)).*hamming(size(T,1));
    % %%% Clustering based on data filtered by the detection filters so that 
    % %%% the spectrum is appropriately normalized 
    % %     Q=[];
    % %     for k = 1:length(hos)
    % %         Q(:,:,k) = fftshift(hos(k).xfilt(z(Trl(:,clusterpeaks)).*hamming(size(Trl(:,clusterpeaks),1))),1);
    % %     end
    % %      Q = Q(abs(tt)<=tolerance*2,:,:);
    % %      mxQ = squeeze(max(Q));
    % %      Q = zscore(Q);
    % %     Q = reshape(permute(Q,[1 3 2]),size(Q,1)*size(Q,3),size(Q,2));
    % %     
    %     
    %     S = Q*Q';
    %     [u,l] = svd(S);
    %     v = Q'*u*diag(diag(l).^-.5);
    %     vz = zscore(v(:,1:ncomp));
    %      cl = dbscan(vz,.5,50);
    %      
    % %     [cl,V,D] = spectralcluster(zscore(v(:,1:ncomp)),keepc);
    %     % 
    %     % figure,
    %     % scatter3(v(:,1),v(:,2),v(:,3),1,cl)
    %     % axis equal vis3d
    % 
    % %     [crt,~,~,lbl] = crosstab(cl,compno*2.^(0:keepc-1)');
    %     [crt,~,~,lbl] = crosstab(cl,2.^(compno(clusterpeaks)-1));
    %     [~,cl2comp] = max(crt(:,1:keepc)./sum(crt(:,1:keepc),2));
    % %     crt3 = crosstab(cl3,compno*2.^(0:keepc-1)');
    %     crt3 = crosstab(cl3,2.^(compno(clusterpeaks)-1));
    %     [~,cl32comp] = max(crt3(:,1:keepc)./sum(crt3(:,1:keepc),2));
    % 
    %     
    %      clear M Mx
    %      for k = 1:length(crt)-1
    %          M(:,k) = mean(Q(:,cl==k),2);
    %          Mx(:,k) = mean(z(Trl(:,clusterpeaks(cl==k))),2);
    %      end
    %      
    %     %%
    %     fig = figure('units','normalized','position',[0 0 1 1]);
    %     cols = hsv(keepc);
    %     for k = 1:keepc
    %         subplot(3,keepc+1,k)
    %         plot(1000*tt([1 end]),[0 0],'k')
    %         hold on
    %         plot(1000*tt,hos(k).feature,'linewidth',2,'color',cols(k,:))
    %         grid on
    %         ylim(max(abs(ylim))*[-1 1])
    %         title(sprintf('Comp. %i',k))
    %         if k==1
    %              ax = axis;
    %             yh = ylabel(sprintf('HOSD\nComponents'),'fontweight','bold','position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
    %         end
    %         subplot(3,keepc+1,k+keepc+1)
    %         M = mua(Trl(:,clusterpeaks(cl==cl2comp(k))));
    %         rp = randperm(size(M,2));
    %         plot(1000*tt,M(:,rp(1:min(100,length(rp)))),'k')
    %         hold on
    %         plot(1000*tt,mean(M,2),'color',cols(k,:),'linewidth',2)
    %         % title(sprintf('Waveform cluster %i (N=%i)',k,size(M,2)))
    % 
    %         if k==1
    %         ax = axis;
    %         th = title(sprintf('Waveform\nclustering'),'position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
    %         end
    % 
    %         subplot(3,keepc+1,k+keepc*2+2)
    %         M3 = mua(Trl(:,clusterpeaks(cl3==cl32comp(k))));
    %         rp = randperm(size(M3,2));
    %         plot(1000*tt,M3(:,rp(1:min(100,length(rp)))),'k')
    %         hold on
    %         plot(1000*tt,mean(M3,2),'color',cols(k,:),'linewidth',2)
    % 
    %         % title(sprintf('Skewness cluster %i (N=%i)',k,size(M3,2)))
    %          if k==1
    %             ax = axis;
    %             th = title(sprintf('Skewness\nclustering'),'position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
    %          end
    %     end
    %     subplot(3,keepc+1,2*keepc+2)
    %     if size(v,2)<3
    %         v(:,3) = 0;
    %     end
    %     plh = scatter3(zscore(v(:,1)),zscore(v(:,2)),zscore(v(:,3)),1,cols(cl2comp(cl),:));
    %     xlabel PC1;ylabel PC2; zlabel PC3;
    %     axis equal vis3d
    %     title('Waveform clusters')
    %     subplot(3,keepc+1,keepc*3+3)
    %     if size(v3,2)<3
    %         v3(:,3) = 0;
    %     end
    %     plh = scatter3(zscore(v3(:,1)),zscore(v3(:,2)),zscore(v3(:,3)),1,cols(cl32comp(cl3),:));
    %     xlabel PC1;ylabel PC2; zlabel PC3;
    %     axis equal vis3d
    %     title('Skewness clusters')
    % 
    %     subplot(3,keepc+1,keepc+1)
    %     xx=[crt(cl2comp,:);crt3(cl32comp,:)];
    %     [I,J] = ndgrid(1:size(xx,1),1:size(xx,2));
    %     imagesc(xx)
    %     set(gca,'xtick',1:size(cl,2))
    %     xl = arrayfun(@(k)sprintf('%i',k),1:keepc,'uniformoutput',false);
    %     bb=dec2bin(1:2^length(xl)-1)=='1';
    %     bb = bb(:,end:-1:1);
    %     xll = arrayfun(@(k)sprintf('%s + ',xl{bb(k,:)}),1:size(bb,1),'uniformoutput',false);
    %     xll = cellfun(@(x)x(1:end-3),xll,'uniformoutput',false);
    %     xlinds = str2double(lbl(:,2));
    %     ylbl = [arrayfun(@(k)sprintf('wave clust. %i',k),1:size(crt,1),'uniformoutput',false),arrayfun(@(k)sprintf('skew clust. %i',k),1:size(crt,1),'uniformoutput',false)];
    %     set(gca,'ytick',1:2*size(crt,1),'xtick',1:size(crt,2),'xticklabel',xll(xlinds),'xticklabelrotation',0,'yticklabel',ylbl,'XAxisLocation','top')
    %     xlabel('HOSD component')
    %     for k = 1:numel(xx)
    %         text(J(k),I(k),sprintf('%i',xx(k)),'horizontalAlignment','center','fontsize',6);
    %     end
    %     title('Component-cluster correspondence')
    %     ylabel('clusters')
    %     xlabel('HOSD components')
    % %%
        clear out
        out.spike_times = pkts/dat.fs(1);
        out.spike_times_aligned = (pkts+dt)/dat.fs(1);
        out.components_number= compno;
    %     out.waveform_cluster(clusterpeaks) = cl; %%% This is clustering based on PCA applied to aligned waveforms.
    %     out.skewness_cluster(clusterpeaks) = cl3; %%% This is clustering based on peri-spike time skewness within the respec
        out.skewness_cluster= cl3; %%% This is clustering based on peri-spike time skewness within the respec
        out.skewness_cluster_means = M3x;
        out.sampling_rate = dat.fs(1);
        out.HOSD_waveforms = [hos(1:keepc).feature];
        out.HOSD_detection_filters = [hos(1:keepc).filterfun];
        out.hosobj = hosminimal(hos(1:keepc));
        out.source_file = datafile;
        out.source_variable = vname;

        [pth,fn] = fileparts(datafile);
        outd = fullfile(outputdir,fn);
        if ~exist(outd,'dir')
            mkdir(outd)
        end


        outputfile = fullfile(outd,sprintf('HOSD_spike_sort_%s',vname));
        savefig(fig1,[outputfile,'.fig'],'compact'                                                                                                                                                                                                                                                                                                                                                                                                                                              );
        set(fig1,'renderer','painters')                                                                                                                             
        set(fig1,'units','normalized','position',[0 0 1 1]);
        print(fig1,[outputfile,'.eps'],'-r400','-depsc')
        save(outputfile,'-struct','out')
        close all

    end
end
%%
% wb = cellfun(@fftshift,hos(1).freqindx.Bfreqs,'uniformoutput',false);
% figure
% for k = 1:keepc
%     subplot(1,keepc,k)
%     imagesc(wb{:},fftshift(abs(hos(k).bicoh)));
%     if k ==1
%         cax = caxis;
%     else
%         caxis(cax)
%     end
%     axis xy
%     
% end
