function blocks=locateNlx(blocks,varargin)

ignore_empty_files  = true;
persistent hbrl_data_dir

interactive = true;

i=1;

while i < length(varargin)

    switch lower(varargin{i})
    
         case 'interactive' 
            interactive = varargin{i+1};
             i = i+1; 
        case 'hbrl_dir'
            hbrl_data_dir= varargin{i+1};
             i = i+1; 
       

    end
    i = i+1;
end

if isempty(hbrl_data_dir) && interactive
    
    hbrl_data_dir=uigetdir('C:\','Please locate the mapped network drive for HBRL_Data.');
    if hbrl_data_dir ==0
        hbrl_data_dir='';
        return
    end
end
   
d = dir(fullfile(hbrl_data_dir,'Data*'));

rd = regexp({d.name},'_(\d*)-(\d*)','tokens','once');
isce = @(x) cellfun(@isempty,x); 
d(isce(rd))=[];
rd =cat(1,rd{:}); 
rd(isce(rd))={'Inf'};
rd = cellfun(@str2double,rd);
 

for k = 1:length(blocks)
   
    subjid = deblank(blocks(k).subject_id);
    
    sdir = d(str2double(subjid)>=rd(:,1) &  str2double(subjid)<=rd(:,2)).name;
    ssdir = dir(fullfile(hbrl_data_dir,sdir,['Subject.',subjid,'*']));    
    ssdir=ssdir.name;
    data_type = regexp(sdir,'Data\d+_([^_]+)','tokens','once');

    basepath =sprintf('%s%s%s%s%s%sData%sNeuraLynx-TDT', hbrl_data_dir,filesep,sdir,filesep,ssdir,filesep,filesep);
    if ~exist(basepath)
        basepath =sprintf('%s%s%s%s%s%sData%sNLX-TDT', hbrl_data_dir,filesep,sdir,filesep,ssdir,filesep,filesep);
        
    end
    dds = dir(basepath);
    dds(strcmp({dds.name},'..'))=[];
%     dds(strcmp({dds.name},'.'))=[];
    for kk = 1:length(dds)
       ddd = dir(fullfile(basepath,dds(kk).name));       
       re = ~cellfun(@isempty,regexp({ddd.name},blocks(k).block));
       if any(re)
           break
       end
    end
    if isempty(dds) || ~any(re)
        error('Directory not found');
    end
    subdir = fullfile(basepath,dds(kk).name,ddd(re).name);
    fn = dir(subdir);
    fn(strcmp({fn.name},'..'))=[];
    fndir = fn([fn.isdir]);
    lfp = {};
    evnt = {};
    inpt = {};
    blkfiles = struct('path','','lfp',[],'evnt',[],'inpt',[],'block','','blkstruct',[]);
    blkfiles = blkfiles([]);
    blocks(k).blkfiles=blkfiles;
    for kk = 1:length(fndir) 
        %%
        ssdir=fullfile(subdir,fndir(kk).name);
        dlfp  = dir(fullfile(ssdir,'LFP*'));
        if ignore_empty_files
           dlfp([dlfp.bytes]<=16384)=[]; 
        end
        if isempty(dlfp)
            continue
        end
        re = regexp({dlfp.name},'LFPx(\d*)_*(\d*)','tokens','once');
        re=strcat('0',cat(1,re{:}));
        chnums = cellfun(@str2num,re);
        [srt,srti] = sortrows(chnums(:,[2 1]));
        unq = unique(srt(:,1));
        for j = 1:length(unq)
            blkfiles(j).path = ssdir;
            blkfiles(j).lfp = {dlfp(srti(unq(j)==srt(:,1))).name};
            if length(unq)<2 && unq==0    
                blkfiles(j).block = blocks(k).block;
            else
                blkfiles(j).block = sprintf('%s_%0.4i',blocks(k).block,unq(j));
            end
            blkfiles(j).blkstruct = blocks(k);
        end
        
        devnt  = dir(fullfile(ssdir,'Events*'));
        re = regexp({devnt.name},'Events_*(\d*)','tokens','once');
        re=strcat('0',cat(1,re{:}));
        evnums = cellfun(@str2num,re);
        [srt,srti] = sortrows(evnums);
        for j = 1:length(unq)
          blkfiles(j).evnt = {devnt(srti(unq(j)==srt(:,1))).name};
        end
        
        dinpt= dir(fullfile(ssdir,'Inpt*'));  
         re = regexp({dinpt.name},'Inpt(\d*)_*(\d*)','tokens','once');
        re=strcat('0',cat(1,re{:}));
       inptnums = cellfun(@str2num,re);     
        [srt,srti] = sortrows(inptnums(:,[2 1]));
        for j = 1:length(unq)
          blkfiles(j).inpt = {dinpt(srti(unq(j)==srt(:,1))).name};
        end
        blocks(k).blkfiles = [blocks(k).blkfiles,blkfiles];
    end
    
end
    
    
    
        