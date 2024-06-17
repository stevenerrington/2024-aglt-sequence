

allfeats = [];

%%% Here to ensure dependency is recognized
 hos=hosobject(3);
 hosmin = hosminimal(hos);
%%

   ld = load(fullfile(inputdir,inputfiles{1}));
   inputdat = ld.inputdat(jobindex);
   re2 = inputdat.re2;
   sid = inputdat.sid;
   cdats = inputdat.cdats;
   tref = inputdat.tref;
   fsref = 1./mean(diff(tref));
   wref = inputdat.wref;
   blk = inputdat.blk;
   getatlas = inputdat.getatlas;
   block = sprintf('%s-%s',blk{:});
   
   outfile = fullfile(outputdir,sprintf('allfeats%03i.mat',jobindex));

   if isempty(cdats) || exist(outfile,'file')
      return
   end
   
   cdatnums = [cdats.Contact_Number];
   
   clear atl
   tic
   nproc = 0;
   nfiles = length(re2);
   outstruc = struct('allfeats',[],'allfeatsphaserand',[]);
   for kk = 1:length(re2)
       
       for atlk = 1:length(getatlas)
           atl(atlk).atlas = getatlas{atlk};
           atl(atlk).label = 'unassigned';
       end
       try
            ld = load(re2{kk},'hos','fs','chan','dat','block','hosphaserand','randseed');
       catch
           continue
       end
        if ~isfield(ld,'block')
            continue
        end
         try
            ld.hos(1).filterfun;
        catch
            continue
        end
        hos = ld.hos;
        hosparams = struct('order',hos(1).order,'buffersize',hos(1).buffersize, 'lowpass',hos(1).lowpass, 'highpass',hos(1).highpass,'glowpass',hos(1).glowpass,'xlowpass',hos(1).xlowpass,'xhighpass',hos(1).xhighpass,'sampling_rate',hos(1).sampling_rate);
        sid = ld.block.subject;
        cnum = ld.chan.contact;
        cdat = cdats(find(cdatnums==cnum,1));
    %     txt = mw.readPage(sprintf('Subject/%s/Contact_%i',sid,cnum));
%          xfilt = ld.hos.xfilt(ld.dat);
%          sk = skewness(xfilt);
    %     mni = cellfun(@str2double,regexp(txt,'CIT168toMNI\s*\|\s*([-.\d])*\|\s*([-.\d])*\|\s*([-.\d])*\s*}}','tokens','once'));
        if ~isempty(cdat)
            mni = [cdat.At_CIT168toMNI_X,cdat.At_CIT168toMNI_Y,cdat.At_CIT168toMNI_Z];
            if isempty(cdat.DKT_label)
                dkt = 'unassigned';
            else
                dkt = cdat.DKT_label.fulltext;
            end
            
            if isempty(cdat.Destrieux_label)
                dest = 'unassigned';
            else
                dest = cdat.Destrieux_label.fulltext;
            end
            
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
  
        for sufx = {'','phaserand'}
            fldn = ['hos',sufx{1}];
            clear features
            
            x = ld.dat;
            if strcmp(sufx,'phaserand')
                reseed(ld.randseed)
                x(isnan(x))=0;
                x = zscore(real(ifft(exp(2*pi*1i*rand(size(x))).*abs(fft(x)))));
            end

             xfilt = ld.(fldn).xfilt(x);
             sk = skewness(xfilt);

            for ii = 1:length(ld.(fldn))

                feat=ld.(fldn)(ii).feature;
                scale = sqrt(mean(feat.^2));
                feat =  feat./scale;

                pdfilt=ld.(fldn)(ii).filterfun;
                pdfscale = sqrt(mean(pdfilt.^2));
                pdfilt =  pdfilt./pdfscale;

                featfilt = ld.(fldn)(ii).xfilt(ld.(fldn)(ii).feature);
                featmx = max(featfilt); % This is a rough measure of SNR
                featFFTnrm = fftshift(abs(fft(featfilt)));

                featnrm = ifft(fft(feat).*abs(pdfilt));
                featnrmscale = sqrt(mean(featnrm.^2));
                featnrm = featnrm./featnrmscale;

                featFFTnrmscale = sqrt(mean(featFFTnrm.^2));
                featFFTnrm = featFFTnrm./featnrmscale;

                fftinterp = interp1(w,featFFTnrm,wref)';
                fftinterp(wref<ld.(fldn)(ii).highpass | wref>ld.(fldn)(ii).lowpass)=nan;

    %             finterp = interp1(t,feat,tref)';
                finterp = ftresamp(feat,ld.fs,fsref,length(tref));
                finterp(tref<min(t)|tref>max(t))=nan;

                fnrminterp = ftresamp(featnrm,ld.fs,fsref,length(tref));
                fnrminterp(tref<min(t)|tref>max(t))=nan;

    %             pdfinterp = interp1(t,pdfilt,tref)';
                pdfinterp = ftresamp(pdfilt,ld.fs,fsref,length(tref));
                pdfinterp(tref<min(t)|tref>max(t))=nan;

                features(ii) = struct('feature',finterp,'feature_scale',scale,'pdfscale',pdfscale,'ftFFTnorm',fftinterp,'ftFFTnorm_scale',...
                                      featFFTnrmscale,'pdfilt',pdfinterp,'featnorm',fnrminterp,'featnorm_scale',featnrmscale,'featfiltmax',featmx,'chan',ld.chan,'contact',ld.chan.contact,'protocol',ld.block.subprotocol,...
                                      'subject',ld.block.subject,'blockid',ld.block.block,'component',ii,...
                                      'CIT168toMNI',mni,'file',re2{kk},'DKT',dkt,'Destrieux',dest,atlbl{:},'hosparams',hosparams,'fs',fsref,'fsorig',ld.fs,'norig',ld.(fldn)(1).buffersize,'skewness',sk(ii));

            end
            outstruc.(['allfeats',sufx{1}]) = [outstruc.(['allfeats',sufx{1}]) ,features];
        end
        nproc = nproc+1;

        
        fprintf('\nProcessed %i of %i (%0.1f%%) %0.1f min remaining',nproc,nfiles,nproc/nfiles*100,toc*(nfiles/nproc-1)/60);

   end
save(outfile,'-struct','outstruc');
