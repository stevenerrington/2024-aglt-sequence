
function blk_data = parse_block_url(url)

% blk_data = parse_block_url(url)
%
% This function reads data from the Experimental Protocol Database web page
% directly into a matlab structure. 
%
% Usage:  (1) navigate to the page for the block you wish to retrieve in
% a browser. (2) Copy and paste the url from the browser adress bar. (3) run
% parse_block_url. (4) paste url into dialog box. Alternatively, url can be 
% passed as an argument.
%
% blk_data is a struct with the following fields
%     .subject_id  - self explanatory
%     .block     -   block identifier
%     .record_id  -   number in the database for the retrieved record.
%     .memo       -   contents of the memo field.
%     .lozchannels  - struct array for low Z channels, with the following fields
%            .contact - contact number
%            .label   -  contact label 
%            .number  - contact number within the group (eg Frontal grid 1, etc)  
%            .channel - recording channel number. 
%     .hizchannels - as above for high Z channels
%     .tables - struct array with unparsed and semi-parsed html for tables
%         .id  - Table ID
%         .body - Table contents
%         .fields - struct with fields and values same as the entries in the table.
%                       eg  .Left_attenuation: 20 dB
%      .get_loz([search term]) - returns all low Z channels whose label includes the search term. 
%      .get_hiz([search term]) - returns all high Z channels whose label includes the search term. 
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/parse_block_url.m $
% $Revision: 789 $
% $Date: 2017-01-29 20:43:19 -0600 (Sun, 29 Jan 2017) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2011
if nargin < 1 || isempty(url)
    url = inputdlg(sprintf('Please enter the full URL for the block you                \nwish to retrieve.'));
end

if ~iscell(url)
    url = {url};
end
sort_by_channel_number = false;

chstruc = struct('contact',[],'label',[],...
                    'number',[],'channel',[]);
                
blk_data = struct('subject_id',[],...
          'block',[],...
          'tank',[],...
      'record_id',[],...
      'protocol',[],...
      'subprotocol',[],...
           'memo',[],...
            'url',[],...
    'lozchannels',chstruc,...
    'hizchannels',chstruc,...
         'tables',chstruc,...
        'get_loz',[],...
        'get_hiz',[],...
        'identifier','');
        
for b = 1:length(url)



    dat = urlread(url{b});


    decell = @(x)[x{:}];

%     blk_data(b).subject_id = decell(regexp(dat,'tbPatientID=(\w*)','tokens','once'));
    blk_data(b).subject_id = decell(regexp(dat,'SubjectID[^>]*value="(\w*)"','tokens','once'));
    blk_data(b).block = decell(regexp(dat,'tbBname[^>]*value="(.*?)"','tokens','once'));
     blk_data(b).identifier = decell(regexp(dat,'tbFName[^>]*value="(.*?)"','tokens','once'));
    blk_data(b).tank = decell(regexp(dat,'tbDTName[^>]*value="(.*?)"','tokens','once'));
    blk_data(b).record_id = str2num(decell(regexp(dat,'RecordID=(\d*)','tokens','once')));
    blk_data(b).date = decell(regexp(dat,'tbExpDate[^>]*value="(.*?)"','tokens','once'));
    blk_data(b).memo =decell(regexp(dat,'tbMemo.*>(.*?)</textarea','tokens','once'));
    blk_data(b).url = url{b};

    xx = [decell(regexp(dat,'select name="dlProtocol"(.+?)</select','tokens','once')),''];
    pn = decell(regexp(xx,'<option selected[^>]*>([^<]+?)</option>','tokens')); 
    blk_data(b).protocol = decell(pn);

    xx = [decell(regexp(dat,'select name="dlSubprotocol"(.+?)</select','tokens','once')),''];
    sbpn = decell(regexp(xx,'<option selected[^>]*>([^<]+?)</option>','tokens')); 
    if strcmp(lower(sbpn),'type in a value')
       sbpn = [decell(regexp(dat,'<input name="tbEnterSubprotocol"[^>]+?value="(.+?)"','tokens','once')),''];
    end
    if iscell(sbpn)
        sbpn = decell(sbpn);
    end
    blk_data(b).subprotocol = sbpn;
    
    [xx,tables] = regexp(dat,'<table\s*id=\s*"(?<id>[^>]*)"\s*border\s*\>(?<body>.*?)</table','tokens','names');
    



    for i = 1:length(tables)

       tdat = regexp(tables(i).body,'\n\s*(<td>.*?)/><','tokens');
       tdat = [ tdat{:}];
       nums=[];
       for k = 1:length(tdat)
        try
           tdin = regexprep(tdat{k},'-',':');
           name =  lower(regexp(tdin,'>\s*([^>\[\]]+?)<','tokens','once'));
           name = [name{:}];
           if isempty(deblank(name)), continue, end;


           
           value =  regexp(tdin,'value="(.*?)"','tokens','once');
           value = [value{:}];

           contacts =  regexp(tdin,'<td>([\d]+[:\d]*)','tokens','once');
           c2 =  regexp(tdin,'<td>[\d:]*(;\s*[\d]*:[\d]*)','tokens','once');
           contacts = [contacts{:},c2{:}];
           
           numrange = regexp(name,'\d+[-:,]\d+','match');
           numrange = str2num([numrange{:},'']);
           
           namerep = regexprep(name,'[/*\#?,")(:]','_');
           namerep = regexprep(namerep,'[ ]*','_');
           namerep = regexprep(namerep,'['']','');
           namerep = regexprep(namerep,'[_]*micro.*','');
           namerep = regexprep(namerep,'[_]*macro.*','');
           namerep = regexprep(namerep,'_of_','_');
        %   namerep = regexprep(namerep,'_*[\d]*_*$','');
           namerep = regexprep(namerep,'_[_\d]*$','');
           tables(i).fields.(namerep) = value;


           if ~isempty(contacts) && ~isempty( str2num(contacts)) %#ok<*ST2NM>
               contacts = str2num(contacts);
              chanmap = regexp(value,'(\d+):(\d+)[^\d]*->[^\d]*(\d+)\s*:\s*(\d+)','tokens');   
              if isempty(chanmap)
                  chs = regexp(value,'(\d+):(\d+)','tokens');
                  if isempty(chs)
                      chs = str2num(regexp(value,'(\d+[:\d,]*)','match','once'));
                  end
                  cmapto ={};
              else
                  chanmap = cellfun(@(x)reshape(x,length(x)/2,2),chanmap,'uniformoutput',false);
                  chs = cellfun(@(x)x(:,2),chanmap,'uniformoutput',false);
                  cmapto = cellfun(@(x)x(:,1),chanmap,'uniformoutput',false);
                  
     %             cmapto = regexp(value,'(\d+):(\d+)[^\d]*->[^\d]*(\d+)\s*:\s*(\d+)','tokens');
              end
              if isempty(chs) && sort_by_channel_number, continue, end
              
              if iscell(chs)
                  chs = cellfun(@(x) cellfun(@str2num,x),chs,'uniformoutput',false);
                  chs = cellfun(@(x) x(1):sign(diff(x)):x(2),chs,'uniformoutput',false);                 
                  channels = [chs{:}];
              else
                  channels=chs;
              end
              
              if ~isempty(cmapto)
                  cmapto = cellfun(@(x) cellfun(@str2num,x),cmapto,'uniformoutput',false);
                  cmapto = cellfun(@(x) x(1):sign(diff(x)):x(2),cmapto,'uniformoutput',false);
                  cmapto = [cmapto{:}];
                  [~,cmapto] = ismember(cmapto,contacts);
              else
                  cmapto = 1:length(channels);
              end
              
           
              if sort_by_channel_number
                  srtindx = channels;
              else
                  srtindx = contacts;
                  cmapto = 1:length(contacts);
                  channels(end+1:length(srtindx)) = nan;
              end
              
              if length(numrange)<length(srtindx);
                  numrange =1:length(srtindx);
              end

              namerep = repmat({namerep},1,length(srtindx))';
              if (isempty(regexp(name,'micro','once'))|| ~isempty(regexp(name,'dialysis','once'))) &&isempty(regexp(name,'High Impedance','once'))
                blk_data(b).lozchannels(srtindx) = struct('contact',num2cell(contacts(cmapto))','label',namerep(:),...
                    'number',num2cell(numrange)','channel',num2cell(channels(:)));
              else
                   blk_data(b).hizchannels(srtindx) = struct('contact',num2cell(contacts(cmapto)'),'label',namerep(:),...
                       'number',num2cell(numrange)','channel',num2cell(channels(:)));
              end
           end
        catch cerr
            warning(cerr.message)
       end
       end
    end

    blk_data(b).tables = tables;
    
    %%% Add functions for retrieving records using search words
    
    blk_data(b).get_loz= @(str) blk_data(b).lozchannels( ~cellfun( @isempty,...
                regexpi( strcat({blk_data(b).lozchannels.label},'#'),str ) ) );
            
    blk_data(b).get_hiz= @(str) blk_data(b).hizchannels( ~cellfun( @isempty,...
                regexpi( strcat({blk_data(b).hizchannels.label},'#'),str ) ) ); %#ok<*AGROW>
            
            
                
    
    
end


