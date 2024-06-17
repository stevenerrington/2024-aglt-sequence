
function err = request_serv(netfolder, varargin)

% Extract files and deliver them to a network folder for distributed
% processing

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/distcomp/request_serv.m $
% $Revision: 1093 $
% $Date: 2018-10-18 21:56:53 -0500 (Thu, 18 Oct 2018) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 1
    netfolder = '\\sharedstorage02.hpc.uiowa.edu\hbrl2\HBRL_Upload\distrib';
end
timeout = 180; % Time out if none is given in the key file
pausedur = 20;
maxfile = 10;
requestdir = fullfile(netfolder,'request');
loop = true;
tstart = now;
fprintf('\nWaiting for a request on %s ...\ntime elapsed:  %6i min\n',requestdir,0)
hbrl_data_dir = '\\sharedstorage02.hpc.uiowa.edu\hbrl2\HBRL_Data';
fsrep = @(x)regexprep(x,'[\\/]',filesep);
while loop
    d = dir(sprintf('%s%s*.key',requestdir,filesep));
    ndone = 0;
    for i = 1:length(d)
        kfn = fullfile(requestdir,d(i).name);
        if ~exist(kfn,'file')
            continue
        end
        fc = readfile(kfn);
        stat = gettok(fc,'status');
        killage = gettok(fc,'status');
        extractdir = fsrep(fullfile(netfolder,gettok(fc,'extractdir')));
        if isempty(extractdir)
            extractdir = requestdir;
        end
        
        if isempty(killage)
            killage=timeout;
        end
        
        switch stat
            

            case 'open'
                fc = regexprep(fc,'open','taken');
                timestamp = gettok(fc,'timestamp');
                fc = regexprep(fc,timestamp,sprintf('%0.11f',now));
                [a,hostname] = system('hostname');
                if a==0
                    fc = sprintf('%s\nextractserver=%s',fc,strtrim(hostname));
                end
        
                blockfile = fullfile(extractdir,gettok(fc,'blockfile'));
                savefile =gettok(fc,'savefile');
                format = gettok(fc,'format');
                
                try
                    bf = load(blockfile);
                catch
                    pause(120)%In case file hasn't finished transferring
                    bf = load(blockfile);
                end      
                if isfield(bf,'restart') && (ischar(bf.restart) && strcmp(bf.restart,'true') || islogical(bf.restart) && bf.restart)
                    bf.restart = false;
                    save(blockfile,'-struct','bf') 
%                      bafn = sprintf('%s.bat',tempname);
%                      [~,tmpba,~] = fileparts(bafn);
                    rsfn = 'restart_temp';
                    fid = fopen([rsfn,'.m'],'w');
                    [~,tl] = system('tasklist');
                    ttpid = regexp(tl,'\nTTankEng[.]exe\s*(\d*)','tokens','once');
                    pid = strtok(java.lang.management.ManagementFactory.getRuntimeMXBean.getName.char,'@');
                      fprintf(fid,'\nsystem(''taskkill /f /pid %s'')\n',pid);
                      if ~isempty(ttpid) && ~isempty(ttpid{1})
                        fprintf(fid,'\nsystem(''taskkill /f /pid %s'')\n',ttpid{1});                                        
                      end
                      fprintf(fid,'cd(''%s'')\n',cd);
                      fprintf(fid,'pause(5)\n');
                    fprintf(fid,'%s',mfilename);
%                     fprintf(fid,'\nquit');                
                    fclose(fid);
                    
                      system(sprintf('matlab -r %s &\n',strtok(rsfn,'.')));
       
%                     fid2 = fopen(bafn,'w');
%                     fprintf(fid2,'cd %s\n',cd);
%                     fprintf(fid2,'\nmatlab -r %s &\n',strtok(rsfn,'.'));
%          %                     fprintf(fid2,'\nsystem(''taskkill /f /pid %s'')',pid);
%                     
% %                     system(sprintf('matlab -r %s &',strtok(rsfn,'.')))
%                     system(bafn)
                               return
                end
                fid = fopen(kfn,'w');
                fprintf(fid,fc);
                fclose(fid);
                
                if ~isequal(requestdir,extractdir)
                    movefile(kfn,extractdir)
                end
                kfn = fullfile(extractdir,d(i).name);

                [pd,~,err] = pulldata(bf.block,bf.pdargs{:},'format',format,'interactive',false,'hbrl_dir',hbrl_data_dir);                   
                if ~isempty(err)
                    warning('Caught error:\n\t%s',err{1}.message) %#ok<WNTAG>
                     fc = regexprep(fc,'taken','error');
                    timestamp = gettok(fc,'timestamp');
                    fc = regexprep(fc,timestamp,sprintf('%0.11f',now));
                    fc = sprintf('%s\nerrmsg=%s',fc,err{1}.message);
                      fid = fopen(kfn,'w');
                    fprintf(fid,fc);
                    fclose(fid);

                    continue
                end
                
                if isequal(pd.record_id,bf.block)
                   newsavefile =sprintf('%s.%s',strtrim(pd.identifier),format); 
                else
                    newsavefile = savefile;
                end
                    
                movefile(pd.file,fullfile(extractdir,newsavefile))
                
                fc = regexprep(fc,'taken','done');
                timestamp = gettok(fc,'timestamp');
                fc = regexprep(fc,timestamp,sprintf('%0.11f',now));
                fc = regexprep(fc,savefile,newsavefile);
                fid = fopen(kfn,'w');
                fprintf(fid,fc);
                fclose(fid);
                ndone = ndone+1;
%                 return
                tstart = now;
                fprintf('\nWaiting for a request on %s ...\ntime elapsed:  %6i min\n',requestdir,0)


            case {'done','error'} 
                ndone = ndone+1;
            case {'finished','closed'}
                
            otherwise
                if (now-d(i).datenum)*24*60 > killage
                    fn = fullfile(netfolder,gettok(fc,'file'));
                    fprintf('\nDeleting %s because it''s %i minutes old',fn,round((now-d(i).datenum)*24*60))
                    delete(fn)
                    delete(kfn)
                    fprintf('\nWaiting for a request on %s ...\ntime elapsed:  %6i min\n',requestdir,round((now-tstart)*24*60))


                end

        end
    end
    if ndone > maxfile
        while ndone > maxfile
   
            fprintf('\nToo many files in the queue, going to sleep for a while...')
            pause(pausedur)

            d = dir(sprintf('%s%s*.key',requestdir,filesep));
            ndone = 0;
            for i = 1:length(d)
                kfn = fullfile(requestdir,d(i).name);
                fc = readfile(kfn);
                stat = gettok(fc,'status');
                ndone = ndone + strcmp(stat,'done');
            end
        end
        tstart = now;
    end

    if pausedur>0 
        fprintf('\b\b\b\b\b\b\b\b\b\b\b%6i min\n',requestdir,round((now-tstart)*24*60))
        pause(pausedur)
    end
end
 
   

        
        
        
        %%% 

function fc = readfile(fn)

fid = fopen(fn);
fc = fread(fid,'uchar=>char')';
fclose(fid);
        
%%%

function val = gettok(fc,tok)

expr = sprintf('%s',tok);
val = regexp(fc,[expr,'=([^\n]*)'],'tokens','once');
val = [val{:}];      

