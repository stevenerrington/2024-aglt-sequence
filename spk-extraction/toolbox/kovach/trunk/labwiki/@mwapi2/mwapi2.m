
classdef mwapi2 < handle

    % Interface for the labwiki API. 
    
    %C. Kovach 2018
    properties
        url = 'https://saccade.neurosurgery.uiowa.edu/labwiki/';
        api='api.php';
        user = '';
        current_action = '';
        query_params;
        logged_in = false;
        tree=[];
        errmsg='';
        saccade_uid = 'Lab';
        wiki_uid='Lab';
        response;
        current_token='';
    end
    properties (Dependent = true)
       timestamp 
        
    end
    properties (Hidden = true)
        
        cookiejar;
        use_oauth=false;
        botmode = true;
        botname = 'matlab_mwapi_access'
        upw;
        connection;
        
    end
        
    properties (Access = private)
       salt = 'f16abf5a_161a_49e7_b2e8_263949b1eb92';
    end
    
    methods  
        
        function me = mwapi2
            try
              me.login();
            catch err
                me.errmsg = err;
                warning(err.message)
            end
        end
     
        %%%%%
        
        function login(me,recredential)
             
            if nargin < 2
                recredential = false;
            end
          
            persistent  loginid cookies crupw
            
            if recredential
                cookies=[];
                me.cookiejar=[];
            end
           if isempty(me.upw) && isempty(crupw) || recredential
                sacuid = input(sprintf('Enter saccade user [%s]:',me.saccade_uid),'s');
                if ~isempty(sacuid)
                    me.saccade_uid=sacuid;
                end
                sacpw = getpw('Enter saccade password:');

                uupw = me.base64enc(int8(sprintf('%s:%s',me.saccade_uid,sacpw)));
                [cupw,hdr] = encrypt(uupw,'',me.salt);
                crupw =  [me.base64enc(cupw),':::',me.base64enc(hdr)];
                me.upw = crupw;

           else
               me.upw=crupw;
           end
           if isempty(cookies) &&~recredential && ~isempty(me.cookiejar)
               cookies=me.cookiejar;
           end
            if isempty(cookies) 
%                 xml='';
              

                newloginid = input(sprintf('Enter Labwiki login [%s]:',loginid),'s');
                 
                if me.botmode && isempty(newloginid)
                    loginid = [me.wiki_uid,'@',me.botname];

                    loginpw = deblank(me.posthttp('https://saccade.neurosurgery.uiowa.edu/webdrop/LW/mwapipw'));      
                elseif isempty(loginid) || recredential 
                    me.botmode=false;
%                     newloginid = input(sprintf('Enter Labwiki login [%s]:',loginid),'s');
                    if ~isempty(newloginid)
                        loginid = newloginid;
                    end
                    loginpw = getpw('Enter LabWiki Password:');
               end
                %me.query_params=struct('action','query','meta','tokens','type','login','format','json');
                logintoken = me.get_token('login');
                %resp = me.posthttp;
               % rjson = me.fromJSON(resp);
                %logintoken = rjson.query.tokens.logintoken;
               % logintoken = regexprep(logintoken,'+\','%2B%5C');

                me.query_params=struct('action','login','lgname',loginid,'lgpassword',loginpw,'format','json','lgtoken',logintoken);
                resp = me.posthttp;
                me.response = me.fromJSON(resp);
                cookies=me.cookiejar;
            else
                me.cookiejar=cookies;
                me.query_params=struct('action','query','meta','userinfo');
                resp = me.request;
                me.response=resp;
                if resp.query.userinfo.id==0;
                    me.login(true);
                    return
                else
                    me.wiki_uid=resp.query.userinfo.name;
                    me.logged_in=true;
                    return
                end
            end
            if isfield(me.response,'login') && strcmpi(me.response.login.result,'Success')
                me.logged_in=true;
                fprintf('\nSuccess...logged in as %s',me.response.login.lgusername)
                me.wiki_uid = me.response.login.lgusername;
            elseif isfield(me.response,'login')
                
                fprintf('\nLogin %s because %s',me.response.login.result,me.response.login.reason)
            end
           
           
        end
        %%%%%%%
           function resp = askargs(me,conditions,showprops,parameters,offset)
            
        % USE:  
        %    result = me.askargs(search_conditions, show_result_properties, ask_parameters)
        %
        % Make and parse a request via api.php/askargs (see
        % see [me.url,'/api.php?action=help&modules=askargs']
        % for details.
        %
        % INPUT:
        %   
        %   search_conditions: Pairs of search property and search value as
        %                       cell array or as the fields/values of a struct.
        %   show_result_properties:  Properties to show in the result
        %   ask_parameters   :  additional parameters for askargs (see api
        %                        help at the url above)
        % OUTPUT:
        %
        %   result is a struct array with the following fields for each
        %   item returned by the search
        %    .main            :  search result page
        %    ,url             :  page url
        %    .(property_name) :  Value of the property "property_name"
            
            limit = 5e3;
            if isstruct(conditions)
                conditions = struct2cell(conditions);
            end
            if nargin < 4
                parameters = {};
            end
            if nargin < 3
                showprops={};
            end
            if nargin < 5
                offset=0;
            end
            if ~iscell(conditions)
                conds=conditions;
            else
                conditions = reshape(conditions(:),2,length(conditions)/2);
                if strcmpi(conditions{1,1},'Category')
                    conds = sprintf('%s:%s',conditions{:,1});
                else
                    conds = sprintf('%s::%s',conditions{:,1});
                end
                for k = 2:size(conditions,2)
                    if strcmpi(conditions{1,k},'Category')
                        conds = sprintf('%s|%s:%s',conds,conditions{:,k});
                    else
                        conds = sprintf('%s|%s::%s',conds,conditions{:,k});
                    end

                end
            end
            me.query_params=struct('action','askargs','conditions',conds);
            props='';
            if nargin >2 && ~isempty(showprops) %Which properties to print
                showprops = regexprep(showprops,'\s','_');
                if ~iscell(showprops) &&any(showprops=='|')
                    props = showprops;
                elseif ~iscell(showprops)
                    showprops={showprops};
                end
                  props = sprintf('%s',props,showprops{1});
                for k =2 :length(showprops)
                          props = sprintf('%s|%s',props,showprops{k});
                end                            
                me.query_params.printouts=props;
            end
             params=sprintf('limit=%i|offset=%i',limit,offset);
            if nargin >3 && ~isempty(parameters) %Ask parameters
                if ~iscell(parameters) 
                    params = parameters;
                else
                     parameters = reshape(parameters(:),2,length(parameters)/2);
           
                    for k =1 :size(parameters,2)
                              params = sprintf('%s|%s=%s',params,parameters{:,k});
                    end                            
                end
            end
            me.query_params.parameters=params;
             resp = me.request;
             if isfield(resp,'query') && isfield(resp.query,'results')
                 if ~isempty(resp.query.results)
                    result = struct2cell(resp.query.results);
                    result = [result{:}];
                    prouts = [result.printouts];
                    if ~isempty(prouts)
                        prfn = fieldnames(prouts);
                        prcell = [prfn';cellfun(@(x){prouts.(x)},prfn,'uniformoutput',false)'];
                    else
                        prcell={};
                    end
                    rescell = [{'main',{result.fulltext}}',{'url',{result.fullurl}}',prcell];
                    resp = struct(rescell{:});
                 else
                     resp = [];
                 end
             end
             if length(resp)==limit
                resp = [resp,me.askargs(conditions,showprops,parameters,offset+limit)] ;
             end

        end
        %%%%%%%
        function resp=request(me,varargin)
            
            if nargin < 2
                query=me.query_params;
            elseif ischar(varargin{1})
                me.current_action=varargin{1};
                query=varargin(2:end);
            else
                query=varargin{1};
            end
            if iscell(query)
                query = struct(query{:});
            end
            if isfield(query,'action')
                me.current_action=query.action;
            else
                query.action=me.current_action;
            end
            if ~isfield(query,'format');
                query.format='json';
            end
            me.query_params=query;
            me.posthttp;
            resp=me.response;
            
            if isfield(query,'format') && strcmpi(query.format,'json') && ischar(resp)
              resp = me.fromJSON(resp);
            end
        end
           
        %%%%%%%%
        
        function tok=get_token(me,type)
           
            if nargin<2
                type='csrf';
            end
            
            me.query_params=struct('action','query','meta','tokens','type',type,'format','json');
            resp = me.posthttp;
            resp=me.fromJSON(resp);
            tok=resp.query.tokens.([type,'token']);
            tok =  regexprep(tok,'+','%2B');
            tok =  regexprep(tok,'\','%5C');
            me.current_token=tok;
        end
       
        
        function strout = base64enc(me,bytesin) %#ok<INUSL>
            import org.apache.commons.codec.binary.Base64;
            
            if ~isa(bytesin,'int8')
                error('Input to base64 encoding must be an array of 8 bit signed integers')
            end
            strout = char(Base64.encodeBase64(bytesin)');

        end
        function bytesout = base64dec(me,bytesin) %#ok<INUSL>
            import org.apache.commons.codec.binary.Base64;
            if ischar(bytesin)
                bytesin=uint8(bytesin);
            end
            if ~isa(bytesin,'int8')&&~isa(bytesin,'uint8') 
                error('Input to base64 encoding must be an array of 8 bit integers')
            end
            bytesout = Base64.decodeBase64(bytesin)';

        end
                  
        function qstr = create_query(me)
            if isstruct(me.query_params)
                fldn = fieldnames(me.query_params);
                qstr = '';
                for k = 1:length(fldn)
                    qstr=sprintf('%s%s=%s&',qstr,fldn{k},deblank(me.query_params.(fldn{k})));
                end
            elseif iscell(me.query_params)

                qstr = sprintf('%s=%s&',me.query_params{:});
                qstr(end)='';
            else
                qstr=' ';
            end
        end
            
        function rqt = get.timestamp(me) %#ok<MANU>
                   
            tnow = datetime(datevec(now),'timezone','local','format','x');
            utcoff = str2double(char(tnow));


            rqt =  sprintf('%i',round((now - utcoff/24 - datenum('1/1/1970'))*24*3600*1e3));
 
        end
        function strout =sha1(me,strin) 
            import java.security.MessageDigest;
           
            mdg = MessageDigest.getInstance('SHA1');
            strout = me.base64enc(mdg.digest(uint8(strin)));
        end
        
        
        %%%%
        function resp = posthttp(me,url,httprequest)

            if nargin < 2 || isempty(url)
                url=[me.url,me.api];
            end
            
            if isempty(me.cookiejar)
               me.cookiejar = java.net.CookieManager;
            end
            
            java.net.CookieHandler.setDefault(me.cookiejar);
            
            postdata = me.create_query;

            re = regexp(url,'(?<protocol>.*?)://(?<host>.*?)(?<file>/.*)','names');

%             ctx = com.sun.net.ssl.SSLContext.getInstance("SSL");
%             
%             keyStore = java.security.KeyStore.getInstance("JKS");
%             kmf = java.security.KeyManagerFactory.getInstance("JKS");
%            
%             final KeyManagerFactory kmf = KeyManagerFactory.getInstance(KeyManagerFactory
%                     .getDefaultAlgorithm());
%             kmf.init(keyStore, KEY_PASSWORD);
%             final TrustManagerFactory tmf = TrustManagerFactory.getInstance(TrustManagerFactory
%                     .getDefaultAlgorithm());
%             tmf.init(keyStore);
% 
%             tm = trustallcerts;
%             tm = tm.getTM;
%             ctx.init('',tm,java.security.SecureRandom());
           
            
            jurl = java.net.URL(re.protocol,re.host,re.file);
            me.connection=jurl.openConnection;

          
            
            if ~isempty(me.upw)
                re = regexp(me.upw,'(?<cupw>.*):::(?<hdr>.*)','names');
                cupw = me.base64dec(re.cupw)';
                hdr = me.base64dec(re.hdr)';
                
                ctupw = decrypt(cupw,hdr,me.salt);
                me.connection.setRequestProperty('Authorization',sprintf('Basic %s',ctupw));
            end
            

            me.connection.setDoOutput(true);
            me.connection.setRequestMethod('POST');
            me.connection.setDoInput(true);
            me.connection.setRequestProperty( 'Content-Type', 'application/x-www-form-urlencoded'); 
            me.connection.setRequestProperty('charset', 'utf-8');
            me.connection.setRequestProperty( 'Content-Length',num2str(length( postdata )));
            me.connection.setUseCaches(false);
            if nargin >2   && ~isempty(httprequest)
                httprequest = reshape(httprequest(:),2,numel(httprequest)/2);

                  for k = 1:size(httprequest,2)
                    me.connection.setRequestProperty(httprequest{1,k},httprequest{2,k});
                  end   
            end
            
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
            if nargout <1
                me.response = resp;
            end
            me.connection.disconnect;
            
        end
        function out=fromJSON(me,in)
            out=matlab.internal.webservices.fromJSON(in);
        end
        
        function downloadFile(me,filename,url)
           
            import org.apache.commons.io.*
            
            java.net.CookieHandler.setDefault(me.cookiejar);
          
            
            re = regexp(url,'(?<protocol>.*?)://(?<host>.*?)(?<file>/.*)','names');
            jurl = java.net.URL(re.protocol,re.host,re.file);
            me.connection=jurl.openConnection;

            
            [pth,fn,ext] = fileparts(filename);
           if exist(fullfile(pwd,pth),'dir')
                %%% Assume the path is relative in this case
                pth = fullfile(pwd,pth);
            end
            
            if ~isempty(me.upw)
                re2 = regexp(me.upw,'(?<cupw>.*):::(?<hdr>.*)','names');
                cupw = me.base64dec(re2.cupw)';
                hdr = me.base64dec(re2.hdr)';
                
                ctupw = decrypt(cupw,hdr,me.salt);
                me.connection.setRequestProperty('Authorization',sprintf('Basic %s',ctupw));
            end
            
%
%             me.connection.setDoOutput(true);
  %          me.connection.setRequestMethod('POST');
 %           me.connection.setDoInput(true);
%            me.connection.setRequestProperty( 'Content-Type', 'application/x-www-form-urlencoded'); 
   %         me.connection.setRequestProperty('charset', 'utf-8');
    %        me.connection.setRequestProperty( 'Content-Length',num2str(length( postdata )));
%             me.connection.setUseCaches(false);
%             if nargin >2   && ~isempty(httprequest)
%                 httprequest = reshape(httprequest(:),2,numel(httprequest)/2);
% 
%                   for k = 1:size(httprequest,2)
%                     me.connection.setRequestProperty(httprequest{1,k},httprequest{2,k});
%                   end   
%             end
            strin = me.connection.getInputStream();
        
            File = java.io.File(fullfile(pth,filename));
            FileUtils.copyInputStreamToFile(strin,File);
            
            
        end
    end 
    
    methods (Hidden=true)
    
        function create_secret_file(me,filename,secret,pw)
            
            if nargin < 2 || isempty(filename)
                filename = me.keyfile;
            end
            if nargin <3 || isempty(secret)
                secret = inputdlg('Enter secret');
                secret=secret{1};
                if isempty(secret)
                    return
                end
            end
            
            if nargin < 4 || isempty(pw)
                pw1='a';
                pw2='b';
                while ~isequal(pw1,pw2)
                    pw1 = getpw('Enter encryption password:');
                    pw2 = getpw('Enter encryption password again:');
                    if ~isequal(pw1,pw2)
                        fprintf('\nPasswords do not match. Try again')
                    end
                end
                pw=pw1;
            end
            finalpw = me.sha1([me.salt,pw]);
            encrypt(secret,filename,finalpw);
                
        end
        function sec = getsec(me,filename,pw)
            if nargin < 2
                filename = me.keyfile;
            end
            if nargin < 3 || isempty(pw) 
                pw = getpw('Enter decryption password:');            
            end
            finalpw = me.sha1([me.salt,pw]);
            sec = decrypt(filename,'',finalpw);
        end
    end
end
