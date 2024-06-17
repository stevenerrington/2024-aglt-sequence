function  [vabl,FAILED_BLOCK,readerrors] = pulldata(recnums, varargin)

% Extracts data from the data server and block information from the protocol 
% server, saving it locally to matfiles.
%
%
%  blkdat = pulldata
% 
% Opens a browser and requires you to copy and paste search results from 
% the protocol search into a dialog window.
%
% Alternatively, 
%
%  blkdat = pulldata(recnums) 
%
% Takes the record number(s) for the block(s) on the protocol server.
%
%  blkdat = pulldata(recnums) 
%
%  Extracts block informtaion and TDT data for the given record number(s) 
%  directly from the protocol and HBRL data servers. The output blkdat (also
%  saved with the data) is a structure that contains block information, such
%  as subject id, protocol name, block id, tank name as well as contact and 
%  numbers and labels.
%
%  The record numbers are located in the right column of the protocol
%  search page. Extracted data will be put into a folder named according to
%  [subprotocol]_[subject id] saved in mat files named [tankname]_[block_name].mat.
%
%  Data extraction requires the TDT API, hence will only work on a windows PC
%  on which the necessary TDT software is installed.
%
%  To enable the script to get data off the data server, the folder HBRL_Data
%  should be mapped to a directory on the local computer. On running the script,
%  you'll be asked the location of the mapped drive. Data tank files will 
%  be copied over and then deleted locally after extraction.
%
%  pulldata(recnums,format,denoise,searchstring)
%
%  Additional arguments: 
%       format- 'kov' for .kov file format.
%       denoise - 'notch' or apply adaptive notch line noise removal before saving to file.
%               - 'dbt'   Denoising with dbtDenoise, which thresholds the
%                        demodulated band transform
%       searchstring - extract only channels whose label contains the given
%                  string.
%       
%
%  See also GET_PROTOCOL_DATA, SCAN_PROTOCOL_BLOCKS and PARSE_BLOCK_URL


%C Kovach 2012

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/pulldata.m $
% $Revision: 978 $
% $Date: 2018-02-08 16:55:41 -0600 (Thu, 08 Feb 2018) $
% $Author: ckovach $
% ------------------------------------------------

%     decell = @(x) [x{:}];
%     recnums = str2num(decell(inputdlg('Enter record numbers (or range as a:z) to extract')));

if nargin < 1
    recnums = [];
end
% 
format = 'mat';
denoise = 'none';  %'notch','dbt'
getExt = true;
ExtChannels  = [1:8];
searchstr = '';
copy_tank_locally=false;
% channel_type = '';
evntcode='Evnt';
% evntcode='EVNT';
getHimp = true;
getLimp = true;
sf = '';
additional_stores = {};
dbtopts = struct('usedbt',false,'bandwidth',1);
dbtopts.args = {'shoulder',0};
interactive = true;
tankfile = [];
i=1;
block = [];
persistent hbrl_data_dir

while i < length(varargin)

    switch lower(varargin{i})
    
        case 'denoise'
            denoise = varargin{i+1};
            i = i+1;
        case {'format','fmt'}
            format = varargin{i+1};
            i = i+1; 
        case {'searchstring','search'}
            searchstr = varargin{i+1};
            i = i+1; 
        case {'external','ext','get ext'}
            getExt = varargin{i+1};
            if ~islogical(getExt(1)) && getExt(1)~=0
               ExtChannels = getExt;
               getExt = true;
            end
            i = i+1; 
        case {'hi','hiz','get high'}
            getHimp = varargin{i+1};
            i = i+1;
        case {'lo','li','loz','get low'}
            getLimp = varargin{i+1};
             i = i+1; 
        case 'suffix'
            sf = varargin{i+1};
             i = i+1;
        case 'dbt'
            dbtopts = varargin{i+1};
             i = i+1;
        case 'additional stores' 
            additional_stores = varargin{i+1};
             i = i+1; 
         case 'interactive' 
            interactive = varargin{i+1};
             i = i+1; 
        case 'hbrl_dir'
            hbrl_data_dir= varargin{i+1};
             i = i+1; 
        case 'tankfile'
            tankfile= varargin{i+1};
             i = i+1; 
        case 'block'
            block= varargin{i+1};
             i = i+1; 

    end
    i = i+1;
end

% if ~ispc
%    error('Sorry... at present pulldata will only work on a Windows PC with necessary TDT software installed') 
% end



if isempty(hbrl_data_dir) && interactive
    
    hbrl_data_dir=uigetdir('C:\','Please locate the mapped network drive for HBRL_Data.');
    
end

 
d = dir(fullfile(hbrl_data_dir,'Data*'));

rd = regexp({d.name},'_(\d*)-(\d*)','tokens','once');
isce = @(x) cellfun(@isempty,x); 
d(isce(rd))=[];
rd =cat(1,rd{:}); 
rd(isce(rd))={'Inf'};
rd = cellfun(@str2double,rd);
 
fprintf('\nScanning data record(s) from the protocol server.');

if isstruct(recnums)
    bl = recnums;
else
    bl = get_protocol_data(recnums);
end
vabl = bl;


if ~ismember(format,{'kov','mat'})
    error('Format has to be ''kov'' or ''mat''')
end

    
%%
if ~iscell(searchstr)
    searchstr= {searchstr};
end

save_sufx = [searchstr{:}];
% 
% sf = '';
% if nargin >= 5
%     switch deblank(lower(channel_type))
%         case {'hi','hiz'}
%             getLimp =false;
%         sf = '_hiZ';    
%         case {'lo','liz','li'}
%             getHimp = false;
%             sf = '_loZ';
%         case {'ex','ext','inpt'}
%             getHimp = false;
%             getLimp =false;
%             sf = '_exinput';
%         otherwise
%             error('%s is an unrecognized channel type (try ''hi'' or ''lo'')')
%     end
%     
% end
save_sufx = [save_sufx,sf];


FAILED_BLOCK = {};
readerrors = {};

get_tank_from_blockdata = isempty(tankfile);
get_block_from_blockdata = isempty(block);
tic
for i = 1:length(vabl)
   try
    subjid = deblank(vabl(i).subject_id);
    savedirs = strcat(vabl(i).subprotocol,'_',subjid);
    
    
    if get_block_from_blockdata
        block = vabl(i).block;
    end
    if get_tank_from_blockdata
        tank = vabl(i).tank;
   
        sdir = d(str2double(subjid)>=rd(:,1) &  str2double(subjid)<=rd(:,2)).name;
        ssdir = dir(fullfile(hbrl_data_dir,sdir,['Subject.',subjid,'*']));    
        ssdir=ssdir.name;
        data_type = regexp(sdir,'Data\d+_([^_]+)','tokens','once');

        switch data_type{1}(1:2)
            case 'RZ'
                fname = sprintf('%s_%s',tank,block);
                basepath =sprintf('%s%s%s\\%s\\Data\\TDT\\%s\\%s\\%s_%s', hbrl_data_dir,filesep,sdir,ssdir,tank,block,tank,block);
            case 'RX'
                fname = sprintf('%s',tank);
                basepath =sprintf('%s%s%s\\%s\\Data\\TDT\\%s', hbrl_data_dir,filesep,sdir,ssdir,tank);
            case 'NL'
                basepath =sprintf('%s%s%s%s%s%sData%sNeuraLynx-TDT', hbrl_data_dir,filesep,sdir,filesep,ssdir,filesep,filesep);
%                 fname = 
            case 'HP'
                error('HP data extraction has not yet been implemented in this script. Beg, and it might happen.')
            otherwise
                error('Unrecognized data type (should be RZ, RX or HP)')

        end

        fprintf('\n###########################\n\nGetting data for %s, subprotocol ''%s''\n\n',vabl(i).block,vabl(i).subprotocol);
        if copy_tank_locally
            tempdir = fullfile(cd,'temp_data_tank');
            if ~exist(tempdir,'dir')
               mkdir(tempdir)
            end

            tevfile=fullfile(tempdir,[fname,'.tev']);
            if exist(tevfile)
                fprintf('Datatank found locally. Extracting data...\n')
            elseif    exist([basepath,'.tev'],'file')        
                fprintf('\nDatatank located on server, copying files from the server...\n');
                tic
                copyfile([basepath,'.*'],tempdir)
                fprintf('%s\nDone copying, now extracting...\n',toc);
            else
                error('Couldn''t locate datatank at %s.tev',basepath)
            end
        else
              tevfile=fullfile([basepath,'.tev']);
              if ~exist(tevfile,'file')
                error('Couldn''t locate datatank at %s.tev',basepath)
              end
        end
    
    else
        [basepath,tank,~] = fileparts(tankfile);
         tevfile = fullfile(basepath,[tank,'.tev']);
         if copy_tank_locally
             tempdir = fullfile(cd,'temp_data_tank');
%         tevfile = fullfile(tempdir,[tank,'.tev']);
            if ~exist(tempdir,'dir')
               mkdir(tempdir)
            end
            if exist(tevfile)
                fprintf('Datatank found locally. Extracting data...\n')
            elseif isempty(basepath)
            elseif    exist([basepath,'.tev'],'file')        
                fprintf('\nDatatank located on server, copying files from the server...\n');
                tic
                copyfile([basepath,'.*'],tempdir)
                fprintf('%s\nDone copying, now extracting...\n',toc);
            else
                error('Couldn''t locate datatank at %s.tev',basepath)
            end
         elseif ~exist(tevfile,'file')
                error('Couldn''t locate datatank at %s.tev',basepath)
         end             

    end

    
    if iscell(savedirs)
        svdr = savedirs{i};
    else
        svdr = savedirs;
    end
    svdr = regexprep(svdr,'[:@%\s]','_');
    if ~exist(svdr)
        mkdir(svdr);
    end
    
  
    tpath = tevfile;

    
    if ~isempty([searchstr{:}])
        lozch= [];
        hizch = [];
        for ssi = 1:length(searchstr)
            lozch = [lozch,vabl(i).get_loz(searchstr{ssi})];
            hizch = [hizch,vabl(i).get_hiz(searchstr{ssi})];
        end
    else
        ise = @(x)cellfun(@isempty,x)| cellfun(@isnumeric,x);
        isused = @(x)cellfun(@(y)isempty(regexpi([y,''],'not[_ ]used')),x);
        getch = @(x) x( ~ise({x.label}) & isused({x.label}) );
        lozch = getch(vabl(i).lozchannels);
        
        if ~isempty(vabl(i).hizchannels)
            hizch = getch(vabl(i).hizchannels);
        else
            hizch=vabl(i).lozchannels([]);
        end
        vabl(i).lozchannels = lozch;
        vabl(i).hizchannels=hizch;
    end
    
    LimpChannels = [lozch.channel];
    HimpChannels = [hizch.channel];
    
    %Get script contents to save with files
    COM.note='This struct contains file contents for extraction scripts';
    fid = fopen([mfilename,'.m'],'r');
    COM.(mfilename) = fread(fid,'uchar=>char')';
    fclose(fid);
    
    [~,ttank,COM.readtdt] = readtdt(tpath,block,'xxxx',1,[],[],[],interactive); % This just initializes the actx control
    storecodenums = invoke(ttank,'GetEventCodes',0);
    if ~isnan(storecodenums)
        for ii = 1:length(storecodenums);
                storecodes{ii} = invoke(ttank,'CodeToString',storecodenums(ii)); %#ok<*AGROW>
                %if ismember(evnt{1})
        end
    else 
        storecodes={};
    end
     evntcode = storecodes(ismember(lower( storecodes),'evnt'));  %Case varies across subjects
     if ~isempty(evntcode)
            [evnt,ttank] = readtdt(tpath,block,evntcode,1,ttank,[],[],interactive); %#ok<ASGLU>
     else
         evnt = [];
     end
     %         [evnt,ttank] = readtdt(tank,block,'EVNT',0);
        stores =[];
        if ~iscell(additional_stores)
            additional_stores = {additional_stores};
        end
        for k = 1:length(additional_stores)
            [stores.(additional_stores{k}),ttank] = readtdt(tpath,block,additional_stores{k},1,ttank,[],[],interactive);               
        end
        savevars ={'evnt','blkdat','COM'}; %variables to save
        if ~isempty(stores)
            savevars{end+1} = 'stores';
        end
        switch denoise
            case {'notch','old'}
                dn = '_linedn';
            case {'dbt'}
                dn = '_dbtdn';
            otherwise
            dn = '';
        end
        fname = fullfile(svdr,sprintf('%s_%s%s%s.%s',tank,block,save_sufx,dn,format));
        vabl(i).file = fname;
        if strcmp(format,'mat')
            if ~exist(fname,'file')
                blkdat = vabl(i); %#ok<NASGU>
                save(fname,'-v6',savevars{:});
            end
        elseif  strcmp(format,'kov')

            clear ks
            if exist(fname)
                delete(fname)
            end
            ks = kstore(fname);
            if ks.newfile
                blkdat = vabl(i); %#ok<NASGU>
                ks.save(savevars{:});
            end

        else 
            error('unrecognized format')
        end

        if getLimp

             licode = storecodes(ismember(lower( storecodes),'lfpx'));  %Case varies across subjects
            %Low impedence
            for ch = LimpChannels
                [x,ttank] = readtdt(tpath,block,licode,ch,ttank,[],[],interactive);

                switch denoise
                    case {'notch','old'}
                        x.dat = removeLineNoise(x.dat,x.fs(1));
                    case 'dbt'
                        x.dat = dbtDenoise(x.dat,x.fs(1),dbtvar{:});
                    case 'none'
                    otherwise
                        error('Unknown denoising option')
                end
                vname = sprintf('li%i',ch);

                fprintf('\nBlock %s, Channel %s',block,vname);

                eval(sprintf('%s=x;',vname));


            if strcmp(format,'mat')
                save(fname,'-v6','-append',vname);
            elseif  strcmp(format,'kov')                    
                ks.save(vname);
            else
                error('unrecognized format')
            end


                clear(vname);
            end
        end
        %High impedance
        if getHimp 

         hicode = storecodes(ismember(lower( storecodes),'pdes'));  %Case varies across subjects

        for ch = HimpChannels
            [x,ttank] = readtdt(tpath,block, hicode,ch,ttank,[],[],interactive);
%             [x,ttank] = readtdt(tank,block,'PDES',ch,ttank);
             switch denoise
                    case {'notch','old'}
                        x.dat = removeLineNoise(x.dat,x.fs(1));
                    case 'dbt'
                        x.dat = dbtDenoise(x.dat,x.fs(1),dbtvar{:});
                    case 'none'

                 otherwise
                        error('Unknown denoising option')
             end
            if dbtopts.usedbt
                x.dat = dbt(x.dat,x.fs(1),dbtopts.bandwidth,dbtopts.args{:});
            end
            vname = sprintf('hi%i',ch);

            fprintf('\nBlock %s, Channel %s',block,vname);

            eval(sprintf('%s=x;',vname));

            if strcmp(format,'mat')
                save(fname,'-v6','-append',vname);
            elseif  strcmp(format,'kov')                    
                ks.save(vname);
            else
                error('unrecognized format')
            end


            clear(vname);
        end
        end
        %external
        if (~iscell(getExt) && getExt) ||  (iscell(getExt) && getExt{i}{j})
             excode = storecodes(ismember(lower( storecodes),'inpt'));  %Case varies across subjects

            for ch = ExtChannels
                [x,ttank] = readtdt(tpath,block,excode,ch,ttank,[],[],interactive);

                vname = sprintf('ex%i',ch);
                fprintf('\nBlock %s, Channel %s',block,vname);

                eval(sprintf('%s=x;',vname));

                if strcmp(format,'mat')
                    save(fname,'-v6','-append',vname);
                elseif  strcmp(format,'kov')                    
                    ks.save(vname);
                else
                    error('unrecognized format')
                end


                clear(vname);
            end
        end
%             delete(fullfile(tempdir,'*'));
        toc
        if strcmp(format,'kov')
            ks.close; delete(ks); 
        end

    catch err
        FAILED_BLOCK{end+1} = block;
        warning('FAILED TO READ BLOCK %s ... skipping',block)
        readerrors{end+1} = err;
        vabl(i).memo = sprintf('Data extraction FAILED with error:\n%s',err.message);
%             rethrow(err)
   end
        
   
    clear ttank
     
end

disp('== End ==');


