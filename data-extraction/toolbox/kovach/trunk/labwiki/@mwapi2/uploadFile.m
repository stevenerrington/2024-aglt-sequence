function res= uploadFile(me,fromfile,tofile,description,comment,clobber)

% result = lwapi.uploadFile(sourcefile,[targetfile],[description],[comment])
%
%   Upload a file to the wiki.
%
% INPUT:
%
%   sourcefile  - File to be uploaded.
%   targetfile  - File name after uploading. The file will exist on the wiki
%                 under ['File:',targetfile]
%   description - Text description to appear on the file's page
%   comment     - Comment in the revision log.
%
%
% OUTPUT: 
%
%   result  -  Structure containing information about the success or failure 
%               of the transfer and other details.
%

% C. Kovach 2018

if nargin < 6 || isempty(clobber)
    clobber = false; % Ignore warnings and replace any existing version of the file
end
[~,fn,ext] = fileparts(fromfile);

if nargin < 3
    tofile = [fn,ext];
end

tok = me.get_token;

if isempty(me.cookiejar)
   me.cookiejar = java.net.CookieManager;
end

java.net.CookieHandler.setDefault(me.cookiejar);

re = regexp([me.url,me.api],'(?<protocol>.*?)://(?<host>.*?)(?<file>/.*)','names');

boundary=sprintf('--FORM BOUNDARY %s FORM BOUNDARY',me.salt);


jurl = java.net.URL(re.protocol,re.host,re.file);
me.connection=jurl.openConnection;

if ~isempty(me.upw)
    re = regexp(me.upw,'(?<cupw>.*):::(?<hdr>.*)','names');
    cupw = me.base64dec(re.cupw)';
    hdr = me.base64dec(re.hdr)';

    ctupw = decrypt(cupw,hdr,me.salt);
    me.connection.setRequestProperty('Authorization',sprintf('Basic %s',ctupw));
end

fid = fopen(fromfile,'r');
dat = fread(fid,'uchar=>char')';
fclose(fid);


%payload=[sprintf('\n\n--%s',boundary),dat,sprintf('\n--%s\n\n',boundary)];
%me.query_params=struct('action','upload','filename',tofile,'format','json','filesize',num2str(length(dat)),'file',payload,'token',tok);
me.query_params=struct('action','upload','ignorewarnings',num2str(clobber),'filename',tofile,'filesize',num2str(length(dat)),'file',fromfile,'format','json','token',tok);

if nargin > 4 && ~isempty(comment)
    me.query_params.comment = comment;
end
if nargin > 2 && ~isempty(description)
    me.query_params.text = description;
end
%postdata = me.create_query;
% postdata(end)='';

 tok=regexprep(tok,'%5C','\');
 tok=regexprep(tok,'%2B','+');
me.query_params.token=tok;

%  if isequal(postdata(end),'&')
%      postdata(end)='';
%  end
postdata = multipart_form(me.query_params,boundary,dat);
 

reqstr = {'Host','saccade.neurosurgery.uiowa.edu',...
'Connection','keep-alive',...
'Content-Length',num2str(length(dat)),...
'Accept','text/plain, */*; q=0.01',...
'Origin','https://saccade.neurosurgery.uiowa.edu',...
'Content-Type',['multipart/form-data; boundary=',boundary],...
'Accept-Language','en-US,en;q=0.9,de;q=0.8'};
reqstr = reshape(reqstr,2,length(reqstr)/2);
me.connection.setDoOutput(true);
me.connection.setRequestMethod('POST');
me.connection.setDoInput(true);
for k = 1:size(reqstr,2)
    me.connection.setRequestProperty( reqstr{1,k},reqstr{2,k});
end
me.connection.setUseCaches(false);
dout = me.connection.getOutputStream();
dout.write(uint8(postdata),0,length(postdata));
dout.close();

strin = me.connection.getInputStream();

len =me.connection.getHeaderField('Content-Length');
if ~isempty(len)
    len=str2double(len);
    resp = zeros(1,len,'uint8');
    for k = 1:len                
        resp(k)=strin.read; 
    end
    resp = char(resp);
else
    resp=[];
    r=0;
    while r>=0  
          r = strin.read;      
        resp(end+1)=r; 

    end
    resp = char(resp);
end
me.response = resp;
res = me.fromJSON(resp);
me.connection.disconnect;

function output = multipart_form(query,boundary,data)

fldn = fieldnames(query);

output=sprintf('--%s',boundary);
for k = 1:length(fldn)
   
    if strcmpi(fldn{k},'file')
       contentType = 'image/jpeg';
       output=sprintf('%s\nContent-Disposition: form-data; name="%s"; filename="%s"\nContent-Type: %s\n\n%s\n--%s',output,fldn{k},query.(fldn{k}),contentType,data,boundary); 
    else
       output=sprintf('%s\nContent-Disposition: form-data; name="%s"\n\n%s\n--%s',output,fldn{k},query.(fldn{k}),boundary);
    end    
end
