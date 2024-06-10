function [table,url] = ask_labwiki(querystring,varargin)


baseurl = 'https://saccade.neurosurgery.uiowa.edu/labwiki/';
%baseurl = 'https://128.255.85.74/nulabwiki';
defuser = 'Lab';

format = 'json';

if verLessThan('matlab','9')
    [table,url] = ask_labwiki_old(querystring,varargin{:});
    return
end
persistent password  uid


chunk = 5000;
% if nargin >1
%     defuser = user;
%     uid = user;
% end
% if nargin > 2
%     password = passwd;    
% end



% if ~isempty(password)
%     noerr = true;
%     try 
%         webread(baseurl,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password)),' CertificateFilename','');
%     catch
%         noerr = false;
%     end
% %     [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
%      getpasswd = ~noerr;
% else
%     getpasswd = true;
%  
% end
% 
% while getpasswd
%     uid = input(sprintf('\nEnter labwiki username [%s]: ',defuser),'s');
%     if isempty(uid)
%         uid=defuser;
%     end
%     try
%         password = getpw('Enter password: ');
%         noerr=true;
%     catch
%         fprintf('\nEnter password: ')
%         [~,password] = system('read'); 
%         noerr = true;
%         try 
%             webread(baseurl,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password)),' CertificateFilename','');
%         catch
%             noerr = false;
%         end
% 
% %         [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
%         if ~noerr
%             warning('Password, username or other error.');
%         end
%     end
%     getpasswd = ~noerr;
%      
% end

table = [];
if ~iscell(querystring)
    querystring={querystring};
end
fprintf('\nRecords retrieved: ');
prn = fprintf('%5i',0);

for k = 1:length(querystring)
    
     questr = querystring{k};
     resstr = regexp(questr,'\?([^|]*)','tokens');
     resstr = [{'main'},resstr{:}];
     resstr = regexprep(resstr,'[^\w]','_');
     resstr(2,:) = {''};
     str = struct(resstr{:});
%     questr = regexprep(questr,'+','%2B');
%     questr = regexprep(questr,'\s','+');
%     questr = regexprep(questr,'?','%3F');

    
    go = true;
    offset = 0;
    
    while go
        
        switch format
            case 'json'
                json = mwapi(struct('action','ask','query',[querystring{k},sprintf('|limit=%i|offset=%i|link=none|format=list',chunk,offset)],'format','json'),varargin{:});
                if ~isfield(json,'query')
                    error(json.error.query{1})
                end
                if ~iscell(json.query.printrequests)
                    json.query.printrequests = {json.query.printrequests};
                end
                if isempty(json.query.results)
                    warning('No results')
                    return
                end
                json.query.printrequests{1}.key='main'; 
                key = cellfun(@(x)x.key,json.query.printrequests,'uniformoutput',false);
                res = json.query.results;
                fldn = fieldnames(res);
                res2 = cellfun(@(x)res.(x),fldn);
                getpro = ~arrayfun(@(x)isempty(x.printouts),res2);
                if any(getpro)
                     fldn2 = arrayfun(@(x)fieldnames(x.printouts),res2(getpro),'uniformoutput',false);
                    res3=squeeze(struct2cell([res2(getpro).printouts]));
                  %  isstr = cellfun(@isstruct,res3);
                    ftxt = cellfun(@(x)isfield(x,'fulltext'),res3);
                    times = cellfun(@(x)isfield(x,'timestamp'),res3);
                    empt = cellfun(@isempty,res3);
                    res3(ftxt)=cellfun(@(x){x.fulltext},res3(ftxt),'uniformoutput',false);
                    res3(times)=cellfun(@(x){str2double(x.timestamp)},res3(times),'uniformoutput',false);
                    if min(size(res3))==1
                        res3=res3';
                    end
                    issc = cellfun(@(x)isscalar(x)&iscell(x),res3);
                    res3(issc) = [res3{issc}];
                    res3(empt) = {''};                
                    res3 = [{{res2.fulltext}},arrayfun(@(k)res3(k,:),1:length(fldn2{1}),'uniformoutput',false)];
                else
                    res3 = {{res2.fulltext}};
                end
                res4 = [key';res3];
                tab = struct(res4{:});
                           
                
            case 'xml'
                xml = mwapi(struct('action','ask','query',[querystring{k},sprintf('|limit=%i|offset=%i|link=none|format=list',chunk,offset)],'format','xml'),varargin{:});
        
             clear re re2 re3
             re = regexp(xml,'<subject fulltext="([^"]*)"(.*?)</subject','tokens');
             re = cat(1,re{:});
             if isempty(re)
                 break
             end
             re2 = regexp(re(:,2),'<property\s*label="([^"]*)">.*?fulltext="([^"]*).*?</property>','tokens');
             re2 = cellfun(@(x)[regexprep(x{1}(1),'[^\w]','_');x{1}(2)],re2,'uniformoutput',false);
             if isempty([re2{:}])
    %          re2 = regexp(re(:,2),'<property\s*label="([^"]*)"><value[^>]*>([^<]*)</value>','tokens');
                 re2 = cellfun(@(x,y)[cat(1,x{:})',{'main';y}],re2,re(:,1),'uniformoutput',false);
                 re2 = cellfun(@(x)[regexprep(x(1,:),'[^\w]','_');x(2,:)],re2,'uniformoutput',false);
             end
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
      
        end
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
for k = 1:size(vals,2)
    s.(vals{1,k})=vals{2,k};
end