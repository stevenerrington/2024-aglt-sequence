% Subject 522:
% ddir = fullfile('~/hbrl/HBRL_Analyses','Anal17_NLX-TDT_505-','Subject.522A.DBS.STN_TDT','Greenlee_Rohl');
% % File name is 522-002_SPECIALevents.mat. We will look at spike from Inpt_RZ2_chn001-Inpt_RZ2_chn003. The data is already high passed filtered.  Attached are the event time stamps for Speech production. 
% datafile = fullfile(ddir,'522-002_SPECIALevents.mat');
% outputdir = 'output';

% Subject 535:
% ddir = fullfile('~/hbrl/HBRL_Analyses','Anal17_NLX-TDT_505-550','Subject.535A.DBS.STN_TDT','Greenlee_Johari_OR_Timing');
% file name  for STN is 535-001_SPECIALevents.mat , We will look at spikes from Inpt_RZ2_chn001-Inpt_RZ2_chn003. Input 3 is in STN
% datafile = fullfile(ddir,'535-001_SPECIALevents.mat');
% and for Nigra is 535-002_SPECIALevents.mat  we will look at Inpt_RZ2_chn001-Inpt_RZ2_chn003, 
ddir = fullfile('~/hbrl','HBRL_Upload','for_JBerger','JBData','STN','SpikeData')
subdirs = dir(fullfile(ddir,'*-*'));
% datafiles = {fullfile(ddir,'273-002 ','273-002_ch001-003.mat'),fullfile(ddir,'269-001 ','269-001_ch001-003.mat'),fullfile(ddir,'253-002 ','253-002_ch001-003.mat')};
datafiles = arrayfun(@(x)fullfile(x.folder,x.name,sprintf('%s_ch001-003.mat',x.name)),subdirs,'uniformoutput',false);

outputdir = fullfile(ddir,'hosd_sorting');

%%% load the data

for di = length(datafiles):-1:1
    datafile = datafiles{di};

    %%% PARAMETERS FOR HOSD
    windur = .01; % Window duration (sec)
    lowpass = 3000; % lowpass cutoff ( Hz )
    ncomp = 4; % Number of components
    window = 'hamming'; %Window applied after segmentation of the input
    tolerance = .001; % Tolerance for the purpose peak detection (width of smoothing kernel in sec. applied to peak detection output).
    skewness_threshold = .05; %Discard components whose filtered output falls below this threshold. 
    normalize_power = true;
    %%% INITIALIZE THE HOS OBJECT
    hos(ncomp) = hosobject(3); 
    %%
    [pth,fn]= fileparts(datafile);
    blk = regexp(fn,'^\d\d\d-\d\d\d','match','once');

    for chan = 3:-1:1
    %     chan = 1;

        vname = sprintf('Inpt_RZ2_chn%03i',chan);
        dat = load(datafile,vname);
        if isfield(dat,vname)
            dat = dat.(vname);
        else
            ld = load(datafile);
            dat.dat = double(ld.PDesMat(chan,:)');
            dat.fs = ld.srate;
            tmst = fullfile(ddir,'Timestamps',sprintf('%s.mat',regexprep(blk,'-','_')));
            ldts = load(tmst);
            dat.FIDX = ldts.FIDX;
            vname = sprintf('chan00%i',chan);
        end
        hos.initialize(round(windur*dat.fs(1)),dat.fs(1),lowpass,[],[],'window','hamming');

        %%% DENOISING: LINE NOISE AND TRANSIENTS
        figure

        [xdn,~,~,spk] = dbtDenoise(double(dat.dat),dat.fs(1),.25,'make plot',true,'remove spikes',false);
    %     [xdn,~,~,spk] = dbtDenoise(double(dat.dat),dat.fs(1),.25,'make plot',true);

        if ~spk.remove_spikes
             spk.filter  =1;
        end
        %%% Filter in MUA range
        mua = hpfilt(zscore(xdn),[dat.fs(1) 300 ]);
        z = zscore(mua) + 0./(spk.filter>.5);
        if normalize_power
            powsmw=hann(.5*dat.fs);
            powsmw = powsmw./sum(powsmw);
            pow = convn(xdn.^2,powsmw,'same');
            xdnrm = xdn./(sqrt((pow + mean(pow))/2));
            znrm = zscore(mua./pow)+ 0./(spk.filter>.5);
        else
            xdnrm = xdn;
            znrm = z;
        end
        
        %%% APPLY HOSD TO DATA
        [A,B] = hos.get_block(znrm); % A and B are segmented data after and before alignment, respectively.

        %%% Apply the feature identification filter
        xfilt = hos.xfilt(z);
        xfilt(isnan(xfilt))=0;

        %%% function with which to compute skewness
        skewfun = @(x)nanmean(x.^3)./nanmean(x.^2).^(3/2);

        skw = skewfun(xfilt); % Compute skewness on filtered data

    %     keepc = min(find(skw>=skewness_threshold,1,'last')+1,ncomp);
        keepc = min(find(skw>=skewness_threshold,1,'last'),ncomp);

        %%% Filtered and thresholded data
        xthr = hos(1:keepc).xthresh(z);

        %%%%% Get the peaks in the thresholded data %%%%%
        %%% Apply tolerance smoothing
        smxthr = convn(full(xthr),hann(round(2*tolerance*dat(1).fs(1))),'same');
        %%% Here combining over components and smoothing
        [~,pks] = getpeak2(sum(smxthr,2));
        %%% Here separately for each component
        [~,pksep] = getpeak2(smxthr);

        pkts = [];
        % compno = [];
        for k = 1:size(pks,2)
            pkts = [pkts,find(pks(:,k)'==1)];
        end

        pksm = convn(full(pksep==1),ones(round(tolerance*dat.fs(1)),1),'same')>0;

        pkts = unique(pkts,'stable');
    %     compno = pksm(pkts,:)==1;
        [~,compno] = max(smxthr(pkts,:),[],2); 
        compno = compno.*(sum(smxthr(pkts,:),2)~=0);
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
        X = hos.xfilt(z(T).*hamming(size(T,1)));

        Xsk = squeeze(skewfun(X));
        Xsk = Xsk(:,~any(isnan(Xsk)));
        S = Xsk'*Xsk;
        [u3,l3] = svd(S);
        v3 = Xsk*u3*diag(1./sqrt(diag(l3)));
    %     tic, cl3 = spectralcluster(zscore(v3),keepc);toc
    %       tic, cl3 = kmeans(zscore(v3),keepc,'replicates',100);toc
          tic,cl3 = isosplit5_mex(v3');toc

    %     eva3 = evalclusters(v3(:,1:ncomp),'linkage', 'CalinskiHarabasz','klist',2:keepc);
    %     cl3 = eva3.OptimalY;

        %%% Alternative: cluster according to waveforms after common alignment
        %%% This applies HOSD again on the waveforms and the waveforms aligned 
        %%% on the respective component
        hos2 = hosobject(hos(1:keepc));
        hos2(end+1) = hosobject(hos(1));
        for kk = 1:length(hos2)
            hos2(kk).use_adaptive_threshold = false;
        end
        [A,B] = hos2.get_block(z(T));
        kp = ~all(all(isnan(A),2),1);
        A = A(:,:,kp);
        B = B(:,:,kp);

    %     Q = reshape(permute(B,[1 3 2]),size(A,1)*size(A,3),size(A,2));
        %%% This uses the reconstructed feature without realignment
        Q = reshape(permute(B,[1 3 2]),size(B,1)*size(B,3),size(B,2));
    %     Q = z(T);
    %     Q = zscore(Q);

        S = Q*Q';
        [u,l] = svd(S);
        v = Q'*u*diag(diag(l).^-.5);

    %     [cl,V,D] = spectralcluster(zscore(v(:,1:ncomp)),keepc);
    %       tic,[cl,V,D] = kmeans(zscore(v(:,1:ncomp)),keepc,'Replicates',100);toc
    %     eva = evalclusters(v(:,1:2*ncomp),'linkage', 'CalinskiHarabasz','klist',2:keepc);
    %     cl = eva.OptimalY;
        tic,cl = isosplit5_mex(v(:,1:ncomp)');toc
        % 
        % figure,
        % scatter3(v(:,1),v(:,2),v(:,3),1,cl)
        % axis equal vis3d

    %     [crt,~,~,lbl] = crosstab(cl,compno*2.^(0:keepc-1)');
        [crt,~,~,lbl] = crosstab(cl,2.^(compno-1));
        [~,cl2comp] = max(crt(:,1:keepc)./sum(crt(:,1:keepc),2),[],1);
    %     crt3 = crosstab(cl3,compno*2.^(0:keepc-1)');
        crt3 = crosstab(cl3,2.^(compno-1));
        
        [~,cl32comp] = max(crt3(:,1:keepc)./sum(crt3(:,1:keepc),2),[],1);
        keepc0=keepc;
        %%
        close all
        
        fig = figure('Units','normalized','position',[0 0 1 .5],'PaperOrientation','landscape');
        keepc = max([keepc0 max(cl) max(cl3)]);
        cl2comp(end+1:size(crt,1)) = setdiff(1:size(crt,1),cl2comp);
        cl32comp(end+1:size(crt3,1)) = setdiff(1:size(crt3,1),cl32comp);
        
        cols = hsv(keepc);
        for k = 1:keepc
            if k <= keepc0
                subplot(3,keepc+1,k)
                plot(tt([1 end]),[0 0],'k')
                hold on

                plot(tt,hos(k).feature,'linewidth',2,'color',cols(k,:))
                grid on
                ylim(max(abs(ylim))*[-1 1])
                title(sprintf('Comp. %i',k))
                if k==1
                     ax = axis;
                    yh = ylabel(sprintf('HOSD\nComponents'),'fontweight','bold','position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
                end
            end
            if k <= max(cl)
                subplot(3,keepc+1,k+keepc+1)
                M = mua(T(:,cl==cl2comp(k)));
                rp = randperm(size(M,2));
                plot(tt,M(:,rp(1:min(100,length(rp)))),'k')
                hold on
                plot(tt,mean(M,2),'color',cols(k,:),'linewidth',2)
                % title(sprintf('Waveform cluster %i (N=%i)',k,size(M,2)))

                if k==1
                ax = axis;
                th = title(sprintf('Waveform\nclustering'),'position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
                end
            end
            if k<=max(cl3)
            subplot(3,keepc+1,k+keepc*2+2)
            M3 = mua(T(:,cl3==cl32comp(k)));
            rp = randperm(size(M3,2));
            plot(tt,M3(:,rp(1:min(100,length(rp)))),'k')
            hold on
            plot(tt,mean(M3,2),'color',cols(k,:),'linewidth',2)
            end
            % title(sprintf('Skewness cluster %i (N=%i)',k,size(M3,2)))
             if k==1
                ax = axis;
                th = title(sprintf('Skewness\nclustering'),'position',[ax(1)-diff(ax(1:2))/5, mean(ax(3:4))],'rotation',90);
             end
        end
        subplot(3,keepc+1,2*keepc+2)
        if size(v,2)<3
            v(:,3) = 0;
        end
        plh = scatter3(v(:,1),v(:,2),v(:,3),1,cols(cl2comp(cl),:));
        xlabel PC1;ylabel PC2; zlabel PC3;
        axis equal vis3d
        title('Waveform clusters')
        subplot(3,keepc+1,keepc*3+3)
        if size(v3,2)<3
            v3(:,3) = 0;
        end
        plh = scatter3(v3(:,1),v3(:,2),v3(:,3),1,cols(cl32comp(cl3),:));
        xlabel PC1;ylabel PC2; zlabel PC3;
        axis equal vis3d
        title('Skewness clusters')

        subplot(3,keepc+1,keepc+1)
    %     xx=[crt(cl2comp,:);crt3(cl32comp,:)];
        xx=[crt;crt3];
        [I,J] = ndgrid(1:size(xx,1),1:size(xx,2));
        imagesc(xx)
        set(gca,'xtick',1:size(xx,2))
        xl = arrayfun(@(k)sprintf('%i',k),1:keepc,'uniformoutput',false);
        bb=dec2bin(1:2^length(xl)-1)=='1';
        bb = bb(:,end:-1:1);
        xll = arrayfun(@(k)sprintf('%s + ',xl{bb(k,:)}),1:size(bb,1),'uniformoutput',false);
        xll = cellfun(@(x)x(1:end-3),xll,'uniformoutput',false);
        xlinds = str2double(lbl(:,2));
        ylbl = [arrayfun(@(k)sprintf('wave clust. %i',k),1:size(crt,1),'uniformoutput',false),arrayfun(@(k)sprintf('skew clust. %i',k),1:size(crt3,1),'uniformoutput',false)];
        set(gca,'ytick',1:size(xx,1),'xtick',1:size(xx,2),'xticklabel',xll(xlinds(1:size(xx,2))),'xticklabelrotation',0,'yticklabel',ylbl,'XAxisLocation','top')
        xlabel('HOSD component')
        for k = 1:numel(xx)
            text(J(k),I(k),sprintf('%i',xx(k)),'horizontalAlignment','center','fontsize',6);
        end
        title('Component-cluster correspondence')
        ylabel('clusters')
        xlabel('HOSD components')
        if isfield(dat,'FIDX')

            [Tstim,ttstim] = chopper([-1 7],dat.FIDX.time,dat.fs(1));
            [unq,~,unqi] =unique(dat.FIDX.evnt);
            Tstim(Tstim<1)=1;
            Tstim(Tstim>length(z))=length(z);
            ximp = smxthr>0;
            smw = hann(dat.fs*.2);
            smw = smw./sum(smw)*dat.fs(1);
            pksm = convn(pksep==1,smw,'same');
    %         ximpsm = convn(ximp,hann(dat.fs*.2)./sum(hann(dat.fs*.2)),'same');
            fig2 = figure('Units','normalized','position',[0 0 1 1],'PaperOrientation','landscape');
     
            for kk = 1:keepc0
               subplot(2,keepc0,kk)
               imagesc(ttstim,[],ximp(Tstim+(kk-1)*length(z))')
               hold on
               plot(ttstim([1 end]),[1 1]'*find(diff(unqi)>0)','r--')
               axis xy
               ylabel trial
               title(sprintf('HOSD component %i',kk))
                 subplot(2,keepc0,kk+keepc0)
                 hold on
               for stk = 1:length(unq)

                   plh(stk) = plot(ttstim,mean(pksm(Tstim(:,unqi==stk)+(kk-1)*length(z)),2));
               end
               xlim(ttstim([1 end]));
                 ylabel 'Mean firing rate (s^{-1})'
               xlabel 'time (s)'
                 grid on
            end
            legend(plh,strcat({'Event '},num2str((unq))))
            colormap(flipud(gray))
        end
    %%
        out.spike_times = pkts/dat.fs(1);
        out.components_active= compno;
        out.waveform_cluster = cl; %%% This is clustering based on PCA applied to aligned waveforms.
        out.skewness_cluster = cl3; %%% This is clustering based on peri-spike time skewness within the respec
        out.sampling_rate = dat.fs(1);
        out.HOSD_waveforms = [hos(1:keepc).feature];
        out.HOSD_detection_filters = [hos(1:keepc).filterfun];
        out.hosobj = hosminimal(hos(1:keepc));
        out.source_file = datafile;
        out.source_variable = vname;

        fid = fopen([mfilename,'.m'],'r');
        out.COM = fread(fid,'uchar=>char')';
        fclose(fid)
        [pth,fn] = fileparts(datafile);
        outd = fullfile(outputdir,fn);
        if ~exist(outd,'dir')
            mkdir(outd)
        end

        

        outputfile = fullfile(outd,sprintf('HOSD_spike_sort_%s',vname));
        savefig(fig,[outputfile,'.fig']);
        set(fig,'renderer','painters')
        print(fig,[outputfile,'.pdf'],'-r400','-bestfit','-dpdf')
        savefig(fig2,[outputfile,'_raster','.fig']);
        set(fig2,'renderer','painters')
        print(fig2,[outputfile,'_raster','.pdf'],'-r400','-bestfit','-dpdf')
        save(outputfile,'-struct','out')
%         close all

    end
end
