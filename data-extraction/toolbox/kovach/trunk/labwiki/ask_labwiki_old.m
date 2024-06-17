function [table,url] = ask_labwiki(querystring,user,passwd)


baseurl = 'https://saccade.neurosurgery.uiowa.edu/labwiki/';
%baseurl = 'https://128.255.85.74/nulabwiki';
defuser = 'Lab';


persistent password  uid


chunk = 5000;
if nargin >1
    defuser = user;
    uid = user;
end
if nargin > 2
    password = passwd;    
end



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
    try
        password = getpw('Enter password: ');
        noerr=true;
    catch
        fprintf('\nEnter password: ')
        [~,password] = system('read'); 
        [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
        if ~noerr
            warning('Password, username or other error.');
        end
    end
    getpasswd = ~noerr;
     
end

table = [];
if ~iscell(querystring)
    querystring={querystring};
end
fprintf('\nRecords retrieved: ');
prn = fprintf('%5i',0);

for k = 1:length(querystring)
%      questr = regexp(querystring{k},'(.*?)(\|\?.*)','tokens','once');
%     if isempty(questr)
%         questr{1} = querystring{k};
%         questr{2}='';
%     end
%     questr = regexprep(questr,'[','%5B');
%     questr = regexprep(questr,']','%5D');
%     questr = regexprep(questr,':','%3A');
     questr = querystring{k};
     resstr = regexp(questr,'\?([^|]*)','tokens');
     resstr = [{'main'},resstr{:}];
     resstr = regexprep(resstr,'[^\w]','_');
     resstr(2,:) = {''};
     str = struct(resstr{:});
    questr = regexprep(questr,'+','%2B');
    questr = regexprep(questr,'\s','+');
    questr = regexprep(questr,'?','%3F');
%     questr = regexprep(questr,'\|','%0A');

    go = true;
    offset = 0;
    
    while go

        url = sprintf('%sapi.php?action=ask&query=%s|mainlabel=main|limit=%i|offset=%i&format=xml',...
                baseurl,questr,chunk,offset);            
         xml = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
         
         clear re re2 re3
         re = regexp(xml,'<subject fulltext="([^"]*)"(.*?)</subject','tokens');
         re = cat(1,re{:});
         if isempty(re)
             break
         end
         re2 = regexp(re(:,2),'<property\s*label="([^"]*)">.*?fulltext="([^"]*).*?</property>','tokens');
         re2 = cellfun(@(x,y)[cat(1,x{:})',{'main';y}],re2,re(:,1),'uniformoutput',false);
         re2 = cellfun(@(x)[regexprep(x(1,:),'[^\w]','_');x(2,:)],re2,'uniformoutput',false);
         strs = repmat({str},length(re2),1);
         strs=cellfun(@(s,fs)stfld(s,fs),strs,re2,'uniformoutput',false);
%          strs = cellfun(@(x)struct(x{:}),re2,'uniformoutput',false);
         tab = [strs{:}];
      
      
%         url = sprintf('%sindex.php?title=Special:Ask&q=%s%%0A&p=format%%3Dbroadtable%%2Fmainlabel%%3Dmain%%2Flink%%3Dnone%%2Fheaders%%3Dplain&po=%s&offset=%i&limit=500',...
%                 baseurl,questr{:},offset);            
%         htxt = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
        

%         hd = regexp(htxt,'<th\s*class="(.*?)".*?</th>','tokens');
%         hd = [hd{:}];
%         if isempty(hd)
% %             table=[];
%             go = false;
% %             continue
%         else
%             fns = regexprep(hd,'[^\w]','_');
%         %     for  k = 1:length(hd)
%         % 
%         %         tb = regexp(htxt,sprintf('<td[^>]*class="%s"[^>]*>(.*?)</td>',hd{k}),'tokens');
%         %         arg{1,k} = fns{k};
%         %         arg{2,k}= [tb{:}];
%         %     end
%             tb=regexp(htxt,'<td[^>]*>(.*?)</td>','tokens');
%             tb = reshape(tb,length(hd),length(tb)/length(hd))';
%             for k = 1:length(hd)
%                 arg{1,k} =fns{k};
%                 arg{2,k} =[tb{:,k}];
%             end
%             tab=struct(arg{:});
        table=[table,tab];
        go = length(tab)==chunk;
        offset = offset+chunk;
        if length(table)>4096 && length(table)-length(unique({table.main}))>500
            %%% Found bug which causes the response to loop back to beginning if the total is over 5K
            go=false;
            [unq,unqi] = unique({table.main},'stable');
            table = table(unqi);
        end
        fprintf(repmat('\b',1,prn));
        prn =fprintf('%5i',length(table));
%         end
    end
end
if isempty(table)
      fprintf('\nNo matching records...\n')
end

fprintf(repmat('\b',1,prn));
prn =fprintf('%5i\n',length(table));

function s = stfld(s,vals)
for k = 1:length(vals)
    s.(vals{1,k})=vals{2,k};
end