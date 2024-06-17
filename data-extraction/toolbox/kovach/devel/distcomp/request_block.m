
function out  = request_block(blocks,netfolder, varargin)

% Request data from the extraction server.
%
%
% See also PULLDATA

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/distcomp/request_block.m $
% $Revision: 572 $
% $Date: 2015-03-14 13:24:43 -0500 (Sat, 14 Mar 2015) $
% $Author: ckovach $
% ------------------------------------------------


i = 1; 
format = 'mat';
pausedur = 10;
requestdir = fullfile(netfolder,'request');

savedir = cd;
genkey = true;
interrupt_at_first_file = false;
use_uuid = true;
clobber = false;
query_protocol_server=false;
timeout = 240; %number of imnutes to wait  by default
extractdir = fullfile('extract',char(java.util.UUID.randomUUID));
while i < length(varargin)

    switch lower(varargin{i})

        case 'key'
            key = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1;
            genkey = false;
        case {'format','fmt'}
            format = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1; 
        case 'maxfiles'
            maxfile = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1; 
        case 'savedir'
            savedir = fullfile(netfolder,'extract',varargin{i+1});
            varargin(i:i+1) = [];
            i = i-1; 
            
        case 'extractdir' %% specify fixed extractdir if files should not be reextracted
            extractdir =fullfile('extract', varargin{i+1});
            varargin(i:i+1) = [];
            i = i-1; 
            use_uuid = false;
        case {'timeout'}
            timeout = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1; 
        case 'clobber'
            clobber = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1; 
         case 'interrupt'  % Stop as soon as the first block has been retrieved
            interrupt_at_first_file = varargin{i+1};
            varargin(i:i+1) = [];
            i = i-1; 
        otherwise 
            i = i+1;
    end
    i = i+1;
end
pdargs = varargin;
loop = true;
if isnumeric(blocks)  && query_protocol_server
       blocks = get_protocol_data(blocks); 
end

out = struct('file','','keyfile','');

openblocks = true(size(blocks));
oldstatusmsg='';
while loop && any(openblocks)
	statusmsg = '';

     for k =1:length(blocks)


        if genkey
            key = char(java.util.UUID.randomUUID);
        end
        block = blocks(k);
        if ~isnumeric(block)
            fn = strtrim(block.identifier);
        else
            fn = num2str(block);
        end


        keyfile = fullfile(netfolder,extractdir,sprintf('%s.key',fn));
        if ~exist(keyfile,'file')
            keyfile = fullfile(requestdir,sprintf('%s.key',fn));
        end
        if exist(keyfile,'file')
            fc = readfile(keyfile);
            if isempty(fc)
                statusmsg=sprintf('%s\Error reading %s. Skipping\n',keyfile);
                continue
            end
            stat = gettok(fc,'status');
            savefile =fullfile(netfolder,extractdir, gettok(fc,'savefile'));
            timestamp = str2double(gettok(fc,'timestamp'));

            switch stat
                case {'claimed','taken'}
                    statusmsg=sprintf('%s\nRequest for %s was claimed by %s %i minutes ago...\n',statusmsg,gettok(fc,'savefile'),gettok(fc,'extractserver'),round((now-timestamp)*24*60));
                    k = k+1;
                    continue                       
                case {'open'}
                    statusmsg=sprintf('%s\nRequest for %s was sent %i minutes ago...\n',statusmsg,gettok(fc,'savefile'),round((now-timestamp)*24*60));
                    k = k+1;
                    continue                       
                case 'done'
                    statusmsg=sprintf('%s\nRequest for %s has been filled, transferring file to %s...\n',statusmsg,gettok(fc,'savefile'),savedir);
                    fprintf('%s',statusmsg)
                    blockfile = fullfile(netfolder,extractdir,gettok(fc,'blockfile'));

                    if exist(savefile,'file')
                        fc = regexprep(fc,'done','closed');
                        fid = fopen(keyfile,'w');
                        fwrite(fid,fc);
                        fclose(fid);
                        movefile(savefile,savedir)
                        delete(blockfile)
                    else
                        warning('%s wasn''t found. Deleting the key file.',savefile)
                        delete(blockfile)
                        delete(keyfile)
                        continue
                    end
                    out.file = fullfile(savedir,gettok(fc,'savefile'));
                    out.keyfile = keyfile;
                    openblocks(k) = false;
                    if interrupt_at_first_file
%                         loop = false;
                        return
                    else
                        continue
                    end
                case 'closed'
                    if ~clobber
                        openblocks(k) = false;
                        statusmsg=sprintf('%s\nRequest for %s has already been filled. Use ''clobber'' to force\n',statusmsg,gettok(fc,'savefile'));
                    else
                        delete(keyfile)
                    end
                    continue
                case 'error'
                    if ~clobber
                        erm = gettok(fc,'errmsg');            
                        openblocks(k) = false;
                       statusmsg= sprintf('%s\nExtraction of %s failed with error\n\t%s\n',statusmsg,gettok(fc,'savefile'),erm);
                    else
                        delete(keyfile)
                    end
                    continue
                otherwise
                    warning('%s has an unrecognized status of "%s", skipping...',savefile,stat)
                    continue
            end
        else
             statusmsg= sprintf('%s\nSending request for %s.%s...\n',statusmsg,fn,format);
        end
        blockfn = sprintf('blockdat_%s.mat',fn);
        blockfile = fullfile(netfolder,extractdir,blockfn);
        if ~exist(fullfile(netfolder,extractdir),'dir')
            mkdir(fullfile(netfolder,extractdir))
        end
        save(blockfile,'block','pdargs')
        [er1,usr]=system('whoami');
        [er2,host]=system('hostname');
        if er1 ~= 0 
            usr='';
        end
        if er2~=0
            host='';
        end

        kfcont = sprintf('status=open\ntimestamp=%0.11f\nkey=%s\nformat=%s\nsavefile=%s.%s\nextractdir=%s\nblockfile=%s\ntimeout=%i\nrequester=%s@%s',now,key,format,fn,format,extractdir,blockfn,timeout,deblank(usr),deblank(host));

        fid = fopen(keyfile,'w');
        if fid < 0
            error('failed to open file on destination folder')    
        end
        fwrite(fid,kfcont);
        fclose(fid);
     end
     if loop && any(openblocks)
	if ~isequal(statusmsg,oldstatusmsg)
            statusmsg=sprintf('%s\nWaiting...\n',statusmsg);
            bk=sprintf(regexprep(oldstatusmsg,'.','\b'));
       		fprintf('%s%s',bk,statusmsg);
		oldstatusmsg = statusmsg;
	end
	pause(pausedur)
     end
end
if use_uuid
    try
        delete(keyfile)
        delete(blockfile)        
        rmdir(fullfile(netfolder,extractdir))
    catch
        warning('Failed to delete temporary folder.')
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

