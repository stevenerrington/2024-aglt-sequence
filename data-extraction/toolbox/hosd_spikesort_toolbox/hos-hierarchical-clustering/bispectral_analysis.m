datadir = '~/hbrl/HBRL_Upload/distrib/rqserv';

%%% Signal toolbox sometimes doesn't work on the cluster due to license
%%% problems. So relevant files have been copied locally.
addpath('signal_toolbox_functions') 

% local_data_dir = 'output2';
check_reg_anls = true;

% lowpass = 200;
group = 'ALL';
matlabver='2018a';
if ~exist('sufx','var')
    sufx='';
end
if ~exist('cluster','var')
    cluster ='argon';
end
if ~exist('getch_force','var')
    getch_force =[];
end
switch lower(group)
    case 'mtl'
     load_contacts = {'hippocamp','amyg'};
    case 'all'
        load_contacts = true;
end
clear mdl model
if exist('submit_profiles','var')
             
    waitn = arrayfun(@(x)sum(cellfun(@(y)strcmp(y.status,'waiting'),x.xnes)),submit_profiles);
    [~,profilenum] = min(waitn);
    if length(unique(waitn))==1
        profilenum = ceil(rand*length(submit_profiles));
    end
    submit_profile = submit_profiles(profilenum);
    profile = submit_profile.profile;
    nslots = submit_profile.nslots;
    queue = submit_profile.queue;
    cluster = submit_profile.cluster;
end
if ~exist('nslots','var')
    nslots=4;
end
if ~exist('profile','var')
    profile='high_mem';
end
if ~exist('queue','var')
    queue='NEUROSURGERY,CCOM,UI-HM,all.q';
end
if ~exist('reextract_data','var')
    reextract_data=false;
end
if ~exist('datefilt','var')
    datefilt=Inf;
end

if ~exist('skipdone','var')
    skipdone = true;
end

if ~exist('local_data_dir','var')
    local_data_dir = 'output';
end

netfolder = '~/hbrl/HBRL_Upload/distrib/';
use_cluster = true;
clear xnes
xnes={};
local_save_dirs = {};
if ~exist('window_events','var')
   window_events=[];
end
if ~exist('regoptsfile','var')
    regoptsfile={};
end
already_done = [];
for kk = 1:length(blks)
%         fprintf('\b\b\b\b\b\b\b%s',blk.block);
        blkns={};
        blk = blks(kk);
        opts.block = blk;
        lozch = blk.lozchannels;
        sf = @(x)~cellfun(@isempty,regexpi({lozch.label},x));
        decell=@(x)cat(1,x{:});
        sfmulti = @(x)any(decell(cellfun(sf,x,'uniformoutput',false)),1);
        
        if iscell(load_contacts)
            getch_all = lozch(sfmulti(load_contacts));
        else
            getch_all = lozch;
        end
        if isempty(getch_all)
            fprintf('No channels found! skipping')
            continue
        end
        dats={};
        nlx = [];
        savedir = fullfile(local_data_dir,blk.block);
        subblks = dir([savedir,'_*']);
        
         touchfile = fullfile(savedir,'isRunning.txt');
         
         xneopts = opts;
         xneopts.queue=queue;
         xneopts.quotalimit=100;
         xneopts.profile=profile;
         xneopts.nslots = nslots;
         xneopts.rerun_if_aborted = 'yes';

       if exist(touchfile,'file')
            fid = fopen(touchfile,'r');
            txt = fread(fid,'uchar=>char')';
            fclose(fid);
            if strcmp(txt,'submitting') || isempty(txt)
                tfd = dir(touchfile);
                if now-tfd.datenum < 0/24
                    fprintf('\nA process for %s is already running. Skipping...',blk.block)
                    continue
                end
            else
                 if ~exist(txt,'file')
                     tmpxne = xargon('xxx');
                     txt = fullfile(tmpxne.argonpath,'jobs',txt);
                 end
                 if exist(txt,'file')
                    xnld= load(txt);
                    xne = xnld.me;
                    xnes{kk}=xne;
        %             fprintf('\n%s appears to be already running...skipping',blk.block)
                    switch xne.status
                        case {'error','finished'}
                            if opts.do_regression
                                 mdl = model;
                                mdl.event(1).times = regopts.evnt.time(:)';
                                mdl.event(1).evnt= regopts.evnt.evnt(:)';
                                mdl.event(1).Trange = regopts.Trange;
                            else
                                mdl = [];
                            end
                              fidtouch = fopen(touchfile,'w');
                              fprintf(fidtouch,'submitting');
                              fclose(fid);
                            xne = send_to_cluster(xne,mdl,xneopts);
                            if exist('submit_profile','var')

                                submit_profile = submit_profiles(profilenum);
                                xne.profile = submit_profile.profile;
                                xne.nslots = submit_profile.nslots;
                                xne.queue = submit_profile.queue;
                            end


                            if isempty(xne.jobindices)

                                xne.finish();
                                delete(touchfile)
                            else
                                xne.make_bash_script;
                                xne.create_job;
                                xne.submit(true);
                                 fidtouch = fopen(touchfile,'w');
                                   fprintf(fidtouch,fullfile(xne.subpaths.log.local,'xobject.mat'));
                                  fclose(fidtouch);
                            end
                            continue
                        case {'deferred'}
                            fprintf('\n%s was deferred; resubmitting...',blks(kk).block)
                            xne.create_job();
                            xne.submit();
                            continue
                        case {'running','waiting','submitted','open'}
                            fprintf('\n%s already running...skipping',blks(kk).block)
                            continue               
                    end
                end
            end
         
%            pause(rand*20) %There seems to be a multi-second delay in updating the touch file, so adding a random delay to reduce the likelihood of clashes.
         
       end
           
       if exist(touchfile,'file')
          fidtouch = fopen(touchfile,'w');
          fprintf(fidtouch,'submitting');
          fclose(fid);
       end

%         subblks = subblks([subblks.isdir]);
        
        existing_dir = dir(fullfile(local_data_dir,[blk.block,'*']));
        existing_dir = existing_dir([existing_dir.isdir]);
        if isempty(existing_dir)
            mkdir(fullfile(local_data_dir,blk.block))
        end
        existing_dir = dir(fullfile(local_data_dir,[blk.block,'*']));
        existing_dir = existing_dir([existing_dir.isdir]);

        for subdiri = 1:length(existing_dir)
             savedir = fullfile(local_data_dir,existing_dir(subdiri).name);
       
            extracted_fn = dir(fullfile(savedir,'*_ch*.mat'));
            extracted_channels = regexp({extracted_fn.name} ,'ch(\d)*.mat','tokens','once');
            if ~isempty(extracted_channels)
                extracted_channels =[extracted_channels{:}];
                extracted_channels = cellfun(@str2double,extracted_channels);
                unextracted = getch_all(~ismember([getch_all.channel],extracted_channels));
                avail_chan = extracted_channels;
            else
                unextracted = getch_all;
            end
            if ~isempty(unextracted)
                failed = data_extract(blk,unextracted,savedir,subdiri);
            end
        end
%          subblks = dir([savedir,'_*']);
%         subblks = subblks([subblks.isdir]);
%         
%         existing_dir = dir(fullfile(local_data_dir,[blk.block,'*']));
%         existing_dir = existing_dir([existing_dir.isdir]);

       
        
        for subblki = 1:length(existing_dir)
            if exist('submit_profiles','var')
             
                waitn = arrayfun(@(x)sum(cellfun(@(y)strcmp(y.status,'waiting'),x.xnes)),submit_profiles);
                if length(unique(waitn)) ==1
                    profilenum = ceil(rand*length(submit_profiles));
                end
                [~,profilenum] = min(waitn);
                
                submit_profile = submit_profiles(profilenum);
                profile = submit_profile.profile;
                nslots = submit_profile.nslots;
                queue = submit_profile.queue;
                cluster = submit_profile.cluster;
            end
            
            local_save_dir = fullfile(local_data_dir,existing_dir(subblki).name);
            if opts.do_regression
                    d = dir(fullfile(local_save_dir,sufx,regopts.subdir,'*ch*.mat'));
                    donech = arrayfun(@(x)str2double(decell(regexp(x.name,'ch(\d*).*[.]mat','tokens','once'))),d);
                    
                    mdl = model;
                    mdl.event(1).times = regopts.evnt.time(:)';
                    mdl.event(1).evnt= regopts.evnt.evnt(:)';
                    mdl.event(1).Trange = regopts.Trange;
            else
                d = dir(fullfile(local_save_dir,sufx,'*_out.mat'));
                d = [d,dir(fullfile(local_save_dir,sufx,'*_hos.mat'))];
                d = d([d.datenum]>= now - datefilt);
                donech = arrayfun(@(x)str2double(decell(regexp(x.name,'(\d*)_out.mat','tokens','once'))),d);
                donech = [donech;arrayfun(@(x)str2double(decell(regexp(x.name,'(\d*)_hos.mat','tokens','once'))),d)];
                
                mdl=[];
            end
            d2 = dir(fullfile(local_save_dir,'*_ch*.mat'));
            if isempty(d2)
                warning('No data found for %s, skipping...',existing_dir(subblki).name)
                continue
            end
            d2name = regexp({d2.name},'ch(\d*).mat','tokens','once');
            d2name = [d2name{:}];
            if ~isempty(d2name) 

                availch = cellfun(@(x)str2num(x),d2name);
                redoch = find((~skipdone |~ismember([getch_all.channel],donech))&ismember([getch_all.channel],availch) |...
                            ismember([getch_all.channel],getch_force) );

                if isempty(redoch)
                    fprintf('\n%s already done, skipping...',blk.block);
                    local_save_dirs{kk} = local_save_dir;
                    already_done(end+1) = true;
                    continue
                else
                    getch=getch_all(redoch);
                end
            else
                getch = getch_all;
            end
            getfiles = d2(ismember(availch,[getch.channel]));

                
                
             
            
           if ~use_cluster
                 for chi = 1:length(getch)

               
                    bsid = run_one_channel_bsp(dat);

                    bsid.blkdat=blk;

                    bsids(kk,chi)=bsid; 
                end
           else
                
               if ~exist(fullfile(local_save_dir,sufx),'dir')
                 mkdir(fullfile(local_save_dir,sufx))
               end
%                if exist('opts','var')
%                    optionsfile = {fullfile(local_save_dir,sufx,'options.mat')};
%                     save(optionsfile{1},'opts','model');
%                else
%                    optionsfile={};
%                end
%               switch cluster
%                   case 'argon'
%                     xne = xargon('run_one_channel_bsp.m');
%                   case 'neon'
%                     xne = xneon('run_one_channel_bsp.m');
%               end          
%               xneopts = opts;
% %              xneopts.datafiles=[fullfile(local_save_dir,{getfiles.name}),optionsfile,regoptsfile];
              xneopts.nparallel=length(getfiles);
             xneopts.queue=queue;
%              xneopts.quotalimit=75;
             xneopts.profile=profile;
              xneopts.local_save_dir= fullfile(local_save_dir,sufx);
             xneopts.nslots = nslots;
%              xneopts.rerun_if_aborted = 'yes';
%              xneopts.use_matlab_version=matlabver;
              
%               xne.create_job;
               
               
               xne = send_to_cluster(fullfile(local_save_dir,{getfiles.name}),mdl,xneopts);
%               xne.datafiles=[fullfile(local_save_dir,{getfiles.name}),optionsfile,regoptsfile];
%               xne.nparallel=length(getfiles);
%               xne.queue=queue;
%               xne.quotalimit=75;
%               xne.profile=profile;
%               xne.local_save_dir= fullfile(local_save_dir,sufx);
%               xne.nslots = nslots;
%               xne.rerun_if_aborted = 'yes';
%               xne.use_matlab_version=matlabver;
%               
%               xne.create_job;
               if isempty(xne.jobindices) && ~strcmp(xne.status,'none')
                 fidtouch = fopen(touchfile,'w');
                  fprintf(fidtouch,'submitting');
                  fclose(fid);

                   xne.finish();
                   delete(touchfile)
               elseif ~strcmp(xne.status,'none')
                  if exist('java.opts','file')
                      copyfile('java.opts',xne.subpaths.bash.local);
                      copyfile('mlcacerts',xne.subpaths.assets.local);
                  end          
             
                   xne.submit;
%                    frewind(fidtouch)
                  fidtouch = fopen(touchfile,'w');
%                    fprintf(fidtouch,fullfile(xne.subpaths.log.local,'xobject.mat'));
                   fprintf(fidtouch,fullfile(xne.uuid,'log','xobject.mat'));
                  fclose(fidtouch);
%                xne.monitor;
                    xnes{end+1}=xne;
               end
              if exist('submit_profiles','var')
    
                  submit_profile.xnes{end+1} = xne;
                  submit_profiles(profilenum) = submit_profile;  
              end
           end
      end
end
    