function failed = data_extract(blkdat,channels,savedir,subdiri,denoising_pre_pipeline,resamp)

if nargin < 5 
   denoising_pre_pipeline = {};
end

if nargin < 6 
    resamp = [];
end

isdtdt=false;

try

    blkdat = locateNlx(blkdat);
    
    nlxpath = blkdat.blkfiles(subdiri).path;
    
catch
    nlxpath ='';
end
if nargin< 2 || isempty(channels)
    channels = blkdat.lozchannels;
end
rqpath = '~/hbrl/HBRL_Upload/distrib/';
tdtpath = fullfile(rqpath,'rqserv');
if isempty(denoising_pre_pipeline)
    do_svd_dn='';
else
    do_svd_dn = denoising_pre_pipeline{~cellfun(@isempty,regexpi(denoising_pre_pipeline,'svd','once'))};
end
fid = fopen([mfilename,'.m'],'r');
COM = fread(fid,'uchar=>char')';
fclose(fid);
if ~isempty(denoising_pre_pipeline) && strcmpi(denoising_pre_pipeline{1},'dbtDenoise')
    figure
end

    
try

    
% if ~isempty(nlxpath) && isfield(blkdat,'tank') &&  isempty(blkdat.tank)
if ~isempty(nlxpath) && isfield(blkdat.blkfiles,'lfp') && ~isempty(blkdat.blkfiles(1).lfp) %isempty(blkdat.tank)
    
        fprintf('\nExtracting %s channel %3i',blkdat.block,0);
        for k = 1:length(channels)
            
            chan = channels(k);
            fni = ~cellfun(@isempty,regexp(blkdat.blkfiles(1).lfp,sprintf('LFPx%i[._]',chan.channel)));
            if ~any(fni)
                continue
            end
            dat = readncs(fullfile(nlxpath,blkdat.blkfiles(1).lfp(fni)));
            if ~isempty(resamp)
                
                dat.dat = resample(dat.dat,resamp(1),resamp(2));
                dat.fs = dat.fs*resamp(1)/resamp(2);
            end
            dat.chan = chan;
            dat.block = blkdat;
            if ~isempty(denoising_pre_pipeline) && strcmpi(denoising_pre_pipeline{1},'dbtDenoise')

                dat.dat = dbtDenoise(dat.dat,dat.fs,.1,'make plot',true);
                drawnow
                dat.denoised = true;

            else
                dat.denoised = false;
            end
            dat.extraction_script = COM;
            fprintf('\b\b\b%3i',chan.channel)
            if isempty(do_svd_dn)
                origdatafile = fullfile(savedir,sprintf('%s_ch%i.mat',blkdat.block,chan.channel));            
                save(origdatafile,'dat','origdatafile','chan','blkdat')
            else
                origdatafile = fullfile(savedir,sprintf('%s_ch%i_%s.mat',blkdat.block,chan.channel,do_svd_dn));            
                origdatafiles{k} =origdatafile;
                dats(k) = dat;
            end
            
        end
        
        
else
    d = dir(fullfile(tdtpath,[blkdat.block,'.*']));
    if ~isempty(d)
        [~,fn,ext] = fileparts(d(1).name);
    else
        ext = '';
    end
    istdt = false;
    switch ext
        case '.mat'
            ks.load = @(x) load(fullfile(tdtpath,[fn,'.mat']),x);
            fmt = 'li%i';
        case '.kov'
            ks = kstore(fullfile(tdtpath,[blkdat.block,'.kov']),'r');
            fmt = 'li%i';
        otherwise
            tdtdat = readtdtlnx(blkdat.blkfiles);
            fmt = 'LFPx%03i';
            ks.load = @(vn)tdtdat;
            istdt = true;
%             out = request_block(blkdat,rqpath,'format','kov','hiz',false,'restart',true,'savedir',tdtpath);
%             ks = kstore(fullfile(tdtpath,[blkdat.block,'.kov']),'r');
    end
     for k = 1:length(channels)

            chan = channels(k);
            vname = sprintf(fmt,chan.channel);
            try
                ld = ks.load(vname);
                dat = ld.(vname);
            catch 
                continue
            end
            dat.chan = chan;
            dat.block = blkdat;
            
            if ~isempty(resamp)
                dat.dat = resample(double(dat.dat),resamp(1),resamp(2));
                dat.fs = dat.fs*resamp(1)/resamp(2);
            end
            
            if ~isempty(denoising_pre_pipeline) && strcmpi(denoising_pre_pipeline{1},'dbtDenoise')
                
                dat.dat = dbtDenoise(dat.dat,dat.fs,.1,'make plot',true);
                dat.denoised = true;
                
            else
                dat.denoised = false;
            end
            dat.exraction_script = COM;

            fprintf('\b\b\b%3i',chan.channel)
            origdatafile = fullfile(savedir,sprintf('%s_ch%i.mat',blkdat.block,chan.channel));
            if isempty(do_svd_dn)
                save(origdatafile,'dat','origdatafile','chan','blkdat')
                origdatafile = fullfile(savedir,sprintf('%s_ch%i.mat',blkdat.block,chan.channel));            
            else
                origdatafile = fullfile(savedir,sprintf('%s_ch%i_%s.mat',blkdat.block,chan.channel,do_svd_dn));            
                origdatafiles{k} =origdatafile;
                dats(k) = dat;
            end

     end
   
    
end
    if ~isempty(do_svd_dn)
       hpf = min(dat.fs(1)/4 ,200);
       try
           Xhp = hpfilt([dats.dat],[dat.fs(1) hpf]);
       catch err
           warning(err.message)
           for k = 1:length(dats)
            Xhp(:,k) = hpfilt([dats(k).dat],[dat.fs(1) hpf]);
           end
       end
       nrm = diag(std(Xhp).^-1);
       [u,l] = svd(nrm*(Xhp'*Xhp)*nrm);
       ncomp = str2double(regexp(do_svd_dn,'\d*$','match','once'));
       if isnan(ncomp)
           ncomp=1;
       end
       W = nrm*(eye(size(l))-u(:,1:ncomp)*u(:,1:ncomp)')*nrm^-1;
       Xsvd = [dats.dat]*W;
      
       for k = 1:length(dats)
          dat = dats(k);
          dat.dat = Xsvd(:,k);
          dat.referencing = W(:,k);
          chan = dat.chan;
          origdatafile=origdatafiles{k};
          save(origdatafile,'dat','origdatafile','chan','blkdat')
            
       end
      if istdt
          fldn = fieldnames(tdtdat);
          fldn(~cellfun(@isempty,regexp(fldn,'LFPx')))=[];
          for k = 1:length(fldn)
                 dat = tdtdat.(fldn{k});
                 origdatafile = fullfile(savedir,sprintf('%s.mat',fldn{k}));            
                save(origdatafile,'dat','origdatafile','blkdat');
          end
      end  
    end
    failed = false;
catch err
    failed = true;
    rethrow(err)
end