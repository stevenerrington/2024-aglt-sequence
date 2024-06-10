function bl = get_protocol_data(recns,most_recent_only)

% A script for getting data from the protocol server. 
%
% blk_data = get_protocol_data
%
% Opens a browser with the HBRL protocol search page. Copy and paste the
% search result table into the dialog box window. Make sure to copy the whole
% table, excluding the column labels, and click OK. Next, within the list box
% that comes up, select the blocks for which you want to load protocol data,
% or click OK to load all of them.
% 
% blk_data = get_protocol_data(startrecn,endrecn)
%
% Loads block data from the protocol server for all recent records between startrecn and endrecn.
% Record numbers are located to the right on the protocol server search page.
%
% blk_data is a struct with the following fields
%     .subject_id  - Subject ID.
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
% Note that occasionally channel numbers are reported in non-standard ways,
% especially for HiZ channels, which will cause the data not to be parsed
% correctly.
%
% See also PULLDATA PARSE_BLOCK_URL 


% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/TDT/get_protocol_data.m $
% $Revision: 797 $
% $Date: 2017-02-23 22:05:44 -0600 (Thu, 23 Feb 2017) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 13

if nargin < 2  %%% If multiple records for the same block exist, return only the most recent.
    most_recent_only = true;
end
url = 'http://vcruprotocol.neurosurgery.uiowa.edu/append.aspx?RecordID=%i';

if  nargin < 1 || isempty(recns)
	
%     helpdlg(['Find the record number or range you wish to scan and give those as arguments to ',mfilename]) 
    pause(1)
    web -browser cruprotocol.neurosurgery.uiowa.edu/newSearch.aspx
    fig = figure;
    uicontrol('style','text','string','Copy and paste search results from the browser here.','units','normalized','position',[.1 .9 .8 .05],'fontsize',16);
    uicontrol('style','pushbutton','string','OK','units','normalized','position',[.6 .05 .2 .1],'callback',@(a,b)uiresume(fig));
    bx = uicontrol('style','edit','max',2,'units','normalized','position',[.1 .2 .8 .7],'HorizontalAlignment','left','background',[1 1 1]);
    uiwait(fig);
    if ~ishandle(fig)
        return
    end
    str = get(bx,'string');
    close(fig);
    indx = 1;
    for i = 1:size(str,1);
        try
            a = regexp(str(i,:),'[^\t]*','match');
            blk = deblank(a{1});
            prot = deblank(a{3});
            sprot = deblank(a{4});
            dn = datenum(a{5});
            recns(indx) = str2double(a{end}); 
            psp = [prot,': ',sprot];
            psp(end+1:50) = ' ';
            psp(50:end) = [];
            lstr(indx) = {[blk,'   ',psp,' ',datestr(dn,'MM/DD/YY hh:ss')]}; %#ok<*AGROW>
        catch err
            warning('Failed to parse the following line:\n\n\t''%s''\n\nwith error\n\t%s\nSkipping...',str(i,:),err.message)
            indx = indx-1;
        end
        indx = indx+1;
    end
    sel =listdlg('ListString',lstr,'listsize',[500 400],'initialvalue',1:length(lstr),'Name','Select Blocks','PromptString','Choose which blocks to include (default is all):');
    if isempty(sel)
        return
    else
        recns = recns(sel);
    end
    
end


%%
for i = 1:length(recns)
    try
        rr = parse_block_url(sprintf(url,recns(i)));
        fprintf('\nscanning %i',recns(i));
	rr.failed= false;
	recs(i)=rr;
    catch err
        warning('Failed to get %i with error: %s',recns(i),err.message)
        recs(i).identifier = '';
	recs(i).failed=true;
    end
       recs(i).identifier = [recs(i).identifier,''];
end

[recs.failed]
recs= recs(~[recs.failed]);
if most_recent_only
    [~,~,rn] = unique({recs.identifier});
else
    rn = 1:length(recs);
end

ise = @(x)cellfun(@isempty,x)| cellfun(@isnumeric,x);
isused = @(x)cellfun(@(y)isempty(regexpi([y,''],'not[_ ]used')),x);
getch = @(x) x( ~ise({x.label}) & isused({x.label}) );

for i = 1:length(recs)    
    bl(rn(i)) = recs(i);
    bl(rn(i)).lozchannels=getch(bl(rn(i)).lozchannels);
    bl(rn(i)).hizchannels=getch(bl(rn(i)).hizchannels);
end

