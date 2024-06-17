function blkdat = search_blocks_labwiki(protocols)

persistent password  uid
re={};
if nargin > 0
    re = regexp(blocksubject,'(\d*)','tokens');
    re = cat(1,re{:});
end

baseurl = 'https://saccade.neurosurgery.uiowa.edu/labwiki/';

defuser = 'Lab';

if ~isempty(password)
    [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    getpasswd = ~noerr;
else
    getpasswd = true;
 
end

while getpasswd
    uid = input(sprintf('\nEnter labwiki username [%s]: ',defuser),'s');
    if isempty(uid)
        uid=defuser;
    end

    fprintf('\nEnter password: ')
    [~,password] = system('read'); 
    [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    if ~noerr
        warning('Password, username or other error.');
    end
    fprintf('\n')
    getpasswd = ~noerr;
     
end

 rep = @(tx)regexprep(tx,'<.*?>','');
rep2 = @(tx)regexprep(tx,'\s','_');
 
if nargin==0;
    url = sprintf('%sindex.php/Category:Protocol',baseurl);
    htxt = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    prlist = regexp(htxt,'<h2>Subcategories</h2>.*?</table>','match','once');
    prlist = regexp(prlist,'<li>(.*?)</li>','tokens');
    prlist = rep([prlist{:}]);
    pri = listdlg('liststring',prlist,'promptstring','Which protocol...');
    protocols = prlist(pri);
elseif ~iscell(protocol)
    protocols = {protocol};
end
    
sprlist = {};
for k = 1:length(protocols)
    proto = protocols{k};
    url = sprintf('%sindex.php/Category:%s',baseurl,rep2(proto));
    htxt = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    spr = regexp(htxt,'<li>(.*?)</li>','tokens');
    sprlist = [sprlist,spr{:}];
end
sprlist = rep(sprlist);
pri = listdlg('liststring',sprlist,'promptstring','Which subprotocol...');
subprotocols = sprlist(pri);
    
blklist = {};
for k = 1:length(subprotocols)
    sproto = regexprep(subprotocols{k},'\s','+');
    url = sprintf('%sindex.php?title=Special:Ask&offset=0&limit=500&q=%%5B%%5BCategory%%3ARecording+Block%%5D%%5D%%5B%%5BHas+subprotocol%%3A%%3A%s%%5D%%5D&p=format%%3Dbroadtable&po=%%3FHas+protocol%%0A%%3FHas+subject%%0A%%3FHas+event+start%%0A',...
                 baseurl,rep2(sproto));
    htxt = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    bltxt = regexp(htxt,'>Has event start.*?</table>','match','once');
    blks = regexp(bltxt,'<td>(.*?)</td>.*?"Has-protocol">(.*?)</td>.*?"Has-subject">(.*?)</td>.*?"Has-event-start">(.*?)</td>','tokens');
    if isempty(blks)
        warning('No blocks found!')
        continue
    end
    blks =rep(cat(1,blks{:}));
    blks(:,3) = regexprep(blks(:,3),'None','');
    blks(:,5) = subprotocols(k);
    
    blklist = cat(1,blklist,blks(:,[1 2 5 4]));
end
if isempty(blklist)
   return 
end
fmtstr = '%-15s%-20s%-30s%-60s';
head = sprintf(fmtstr,'Block','Protocol','Subprotocol','time');
liststr = [{head};cellfun(@(a,b,c,d) sprintf(fmtstr,a,b,c,d),blklist(:,1),blklist(:,2),blklist(:,3),blklist(:,4),'uniformoutput',false)];
blk = listdlg('liststring',liststr,'promptstring','Choose a Block...','listsize',[500 400]);

blkdat = get_protocol_labwiki(blklist(blk-1,1),uid,password);