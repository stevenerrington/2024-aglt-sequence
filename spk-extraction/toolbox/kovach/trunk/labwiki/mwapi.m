function out = mwapi(request,user,passwd,mwlogin)

% out = mwap(request)
%
%General purpose script for interacting with the wiki.
%
% The input, "request", is a structure with fields that are named 
% according to the parameters described in the api help. The first  should
% be action. 
%
% Example:
%
% The following code appends "Some test to append" to the page called
% "apiTestEdit".
%
%   request = struct('action','edit','title','apiTestEdit',...
%                     'appendtext','Some text to append');  
% 
%   out =  mwapi(request);
%
% See https://saccade.neurosurgery.uiowa.edu/labwiki/api.php
%

 

if verLessThan('matlab','9.0')
   args = {request};
   if nargin > 1
       args{end+1}=user;
       if nargin > 2
           args{end+1} = passwd;
           if nargin > 3 
               args{end+1}=mwlogin;
           end
       end
   end
   out = mwapi_old(args{:});
   return
end

baseurl = 'https://saccade.neurosurgery.uiowa.edu/labwiki/api.php';

if nargin > 0 && ~isstruct(request)
   error('First input must be a structure.') 
end
    
defuser = 'Lab';
if nargin < 4 || isempty(mwlogin)
    mwlogin = false;
end

persistent password token uid mwapipw

if nargin >1
    defuser = user;
    uid = user;
end
if nargin > 2
    password = passwd;    
end



if ~isempty(password)
    noerr = true;
    try 
        webread(baseurl,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename',''));
    catch
        noerr = false;
    end

%     [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
    getpasswd = ~noerr;
else
    getpasswd = true;
 
end


if ~isfield(request,'format')
    request.format = 'json';
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
        noerr = true;
        try 
            webread(baseurl,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename',''));
        catch
            noerr = false;
        end

%         [~,noerr] = urlread(baseurl,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
        if ~noerr
            warning('Password, username or other error.');
        end
    end
    getpasswd = ~noerr;
     
 
end
if isempty(mwapipw) && mwlogin
    
%         mwapipw = urlread('https://saccade.neurosurgery.uiowa.edu/webdrop/LW/mwapipw','Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
          login = struct('action','login','lgname','Lab@matlab_mwapi_access','format','json');
          out = mwapi(login,uid,password,false);
          logintoken = regexp(out,'.*"token":\s*"([^"]*)".*','tokens','once');
            url = sprintf('%s?action=login&lgname=Lab@matlab_mwapi_access&',baseurl,mwapipw);
          login.lgtoken = logintoken;
            res = webread(url,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename','','RequestMethod','post'),...
                         'lgname','Lab@matlab_mwapi_access','lgpassword',mwapipw,'lgtoken',logintoken{1});

%            [res,err] = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password),...
%                          'post',{'lgname','Lab@matlab_mwapi_access','lgpassword',mwapipw,'lgtoken',logintoken{1}});
%        
end

if isempty(token)
    
   url = sprintf('%s?action=query&meta=tokens&format=xml',baseurl);
   xml = webread(url,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename',''));
%    xml = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password));
   token = regexp(xml,'csrftoken="([^"]*)"','tokens','once');
   token = token{1};
end

fldns = setdiff(fieldnames(request),'action','stable');

reqstring = sprintf('action=%s',request.action);

for k = 1:length(fldns)
    req=regexprep(request.(fldns{k}),'=','%3D');
    req=regexprep(req,'+','%2B');
    req=regexprep(req,'?','%3F');
    reqstring = sprintf('%s&%s=%s',reqstring,fldns{k},req);
end

% reqstring = sprintf('%s&token=%s',reqstring,token);

url = sprintf('%s?%s',baseurl,reqstring);
url = regexprep(url,'\n','%0A');
url = regexprep(url,'\\n','%0A');
url = regexprep(url,'\s','%20');
url = regexprep(url,'#','%23');

if isfield(request,'post') && ~isempty(request.post)
   fldns = fieldnames(request.post);
   for k = 1:length(fldns)
       fldns{k,2} = request.post.(fldns{k,1});
   end
   postdata = fldns;
else
    postdata = {'token',token};
end
 out = webwrite(url,'token',token , weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename','','RequestMethod','post'));
% out = webread(url,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename','','contentType','json'));
%  out = webread(url,weboptions('Timeout',120,'UserName',deblank(uid), 'Password',deblank(password),'CertificateFilename','','contentType','text'));
% out = urlread(url,'Authentication','Basic','Username',deblank(uid),'Password',deblank(password),'post',{'token',token});
