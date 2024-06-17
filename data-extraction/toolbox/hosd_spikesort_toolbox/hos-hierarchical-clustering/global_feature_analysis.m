
[res,out] = system('find ~/lss/bispectral_analysis/ -maxdepth 1 -regex .*[0-9]*-[0-9]*');

re = regexp(out,'[^\n]*','match');

[res,out] = system(sprintf('find ~/lss/bispectral_analysis/  -maxdepth 3 -regex .*hos.*[.]mat'));
nfiles = length(regexp(out,'[^\n]*'));

getatlas = {'BN246','KN_Anatomical_region'};
mw = mwapi2;
 wref = 0:.25:150;
 tref = -2:.002:2; 
%  url = 'https://saccade.neurosurgery.uiowa.edu/labwiki/index.php/Subject/';
allfeats = [];
nproc = 0;
tic

%%
for k = 1:length(re)
    
    [res,out] = system(sprintf('find %s -maxdepth 3 -regex .*hos.*[.]mat',re{k}));
  
    re2 = regexp(out,'[^\n]*','match');
%     re3 = regexp(re2,'/([^/]*)/[^/]*$','tokens','once');
%     re3  = [re3{:}];
%     [unq,unqi] = unique(re3);
    blk = regexp(re{k},'(\d\d\d)-(\d\d\d)','tokens','once');
   sid = blk{1};
   cdats = mw.askargs({'Subject',sid,'Category','Electrode Contact'},[{'Contact Number','At CIT168toMNI X','At CIT168toMNI Y','At CIT168toMNI Z','DKT label','Destrieux label'},strcat(getatlas,' label')],{'link','none'});
   if isempty(cdats)
       continue
   end
   
   cdatnums = [cdats.Contact_Number];
   
   clear atl
   for kk = 1:length(re2)
       
       for atlk = 1:length(getatlas)
           atl(atlk).atlas = getatlas{atlk};
           atl(atlk).label = 'unassigned';
       end
       try
            ld = load(re2{kk},'hos','fs','chan','dat','block');
       catch
           continue
       end
        if ~isfield(ld,'block')
            continue
        end
        hos = ld.hos;
        hosparams = struct('order',hos(1).order,'buffersize',hos(1).buffersize, 'lowpass',hos(1).lowpass, 'highpass',hos(1).highpass,'glowpass',hos(1).glowpass,'xlowpass',hos(1).xlowpass,'xhighpass',hos(1).xhighpass,'sampling_rate',hos(1).sampling_rate);
        sid = ld.block.subject;
        cnum = ld.chan.contact;
        cdat = cdats(cdatnums==cnum);
    %     txt = mw.readPage(sprintf('Subject/%s/Contact_%i',sid,cnum));
    %     xfilt = ld.hos.xfilt(ld.dat);
    %     sk = skewness(xfilt);
    %     mni = cellfun(@str2double,regexp(txt,'CIT168toMNI\s*\|\s*([-.\d])*\|\s*([-.\d])*\|\s*([-.\d])*\s*}}','tokens','once'));
        if ~isempty(cdat)
            mni = [cdat.At_CIT168toMNI_X,cdat.At_CIT168toMNI_Y,cdat.At_CIT168toMNI_Z];
            dkt = cdat.DKT_label.fulltext;
            dest = cdat.Destrieux_label.fulltext;
            
            for atlk = 1:length(atl)
                fld = sprintf('%s_label',atl(atlk).atlas);

                if isfield(cdat.(fld),'fulltext')
                    atl(atlk).label = cdat.(fld).fulltext;
%                     bn246 = cdat.BN246_label.fulltext;
%                 else
%                     bn246 = 'unassigned';
                end
            end
        else
%         else
                dkt = 'unassigned';
                dest = 'unassigned';
    %             bn246 = 'unassigned';
                mni = [nan nan nan];
        end 
        atlbl={atl.atlas;atl.label};
        
        if isempty(mni)
            mni = [nan nan nan];
        end
        w = fftshift(ld.hos(1).sampt)/ld.hos(1).buffersize*ld.fs;
        t = fftshift(ld.hos(1).sampt)./ld.fs;
        clear features
        for ii = 1:length(ld.hos)

            feat=ld.hos(ii).feature;
            scale = sqrt(mean(feat.^2));
            feat =  feat./scale;

            pdfilt=ld.hos(ii).filterfun;
            pdfscale = sqrt(mean(pdfilt.^2));
            pdfilt =  pdfilt./pdfscale;

            featfilt = ld.hos(ii).xfilt(ld.hos(ii).feature);
            featmx = max(featfilt); % This is a rough measure of SNR
            featnrm = fftshift(abs(fft(featfilt)));

            featnrmscale = sqrt(mean(featnrm.^2));
            featnrm = featnrm./featnrmscale;

            fftinterp = interp1(w,featnrm,wref)';
            fftinterp(wref<ld.hos(ii).highpass | wref>ld.hos(ii).lowpass)=nan;

            finterp = interp1(t,feat,tref)';
            finterp(tref<min(t)|tref>max(t))=nan;

            pdfinterp = interp1(t,pdfilt,tref)';
            pdfinterp(tref<min(t)|tref>max(t))=nan;
            
            features(ii) = struct('feature',finterp,'feature_scale',scale,'pdfscale',pdfscale,'ftnorm',fftinterp,'ftnorm_scale',...
                                  featnrmscale,'featfiltmax',featmx,'chan',ld.chan,'contact',ld.chan.contact,'protocol',ld.block.subprotocol,...
                                  'subject',ld.block.subject,'blockid',ld.block.block,'component',ii,...
                                  'CIT168toMNI',mni,'file',re2{kk},'DKT',dkt,'Destrieux',dest,atlbl{:},'hosparams',hosparams);

        end
        allfeats = [allfeats,features];
        nproc = nproc+1;
        fprintf('\nProcessed %i of %i (%0.1f%%) %0.1f hrs remaining',nproc,nfiles,nproc/nfiles*100,toc*(nfiles/nproc-1)/3600);

      end
end
save allfeats allfeats
%% Compute pairwise bispectral correlation using the efficient technique of McGlaughlin 1968
af1 = allfeats([allfeats.component]==1);

mx = [af1.featfiltmax];
[srt,srti] = sort(mx,'descend');
afsrt = af1(srti);

compno = [afsrt.component];
F = [afsrt.feature];
F(isnan(F)) = 0;
PDF = [afsrt.pdfilt];
PDF(isnan(F)) = 0;

ftF = fft(F);
ftPDF = fft(PDF);
% FF = ifft(abs(ftF).^2);
FF = ifft(ftF.*ftPDF);
s3 = sum(FF.^3);
s4 = sum(FF.^4);
K3 = zeros(size(F,2));
K4 = zeros(size(F,2));
for fk = 1:size(F,2)
%     fF = ifft(ftF.*conj(ftF(:,fk)));
    fF = ifft(ftF.*ftPDF(:,fk)); % This implicitly normalizes by SNR.
    K3(fk,:) = sum(fF.^3)./sqrt(s3.*s3(fk));
    K4(fk,:) = sum(fF.^4)./sqrt(s4.*s4(fk));
    fk
end
K=K3;
KK = K(1:fk,:)*K(1:fk,:)';
[u,l] = svd(KK);
v = diag(sqrt(diag(l)).^-1)*u'*K(1:fk,:);

L = linkage(v(1:10,:)','weighted','correlation');
c = cluster(L,'cutoff',.6);

Fft = [afsrt.ftnorm];
Fft(isnan(Fft))=0;
Cfrq = corr(Fft');
Cfrq(isnan(Cfrq))=0;

[ufrq,lfrq] = svd(Cfrq);
vfrq = diag(diag(lfrq).^-.5)*ufrq'*zscore(Fft);

%% SOM classificaton
% x = vfrq(1:20,:);
x = v(1:10,:);

net = selforgmap([4 4]);

% Train the Network
[net,tr] = train(net,x);
y = net(x);
c = (1:size(y,1))*y;
MC = Fft*y'*(y*y')^-1;
MCx = x*y'*(y*y')^-1;
[srtc,srtic] = sort(c);
%% Weighted Phase slope
f = fft(ifftshift(F,1));
ff = f(2:end/2,:).*conj(f(1:end/2-1,:));
phslope = (w(1:end/2-1)*(w(1:end/2-1)'.*diag(abs(ff))))^-1*(w(1:end/2-1))*ff; %Linear slope of group delay vs frequency, weighted by power
%%
x2 = MCx;
net2 = selforgmap([4 4]);
% Train the Network
[net2,tr2] = train(net2,x2);
y2 = net2(x2);
c2 = (1:size(y2,1))*y2;

%% Atlas classifer

atl = {afsrt.BN246};
[unq,~,unqi] = unique(atl);
atlind = zeros(length(unq),size(atl,2));
atlind(unqi'+(0:size(atl,2)-1)*length(unq))=1;

proto = {afsrt.protocol};
[unqpr,~,unqipr] = unique(proto);
protind = zeros(length(unqpr),size(proto,2));
protind(unqipr'+(0:size(proto,2)-1)*length(unqpr))=1;

blockatl = [atlind;protind];

x = vfrq(1:10,:);



