
classdef mwapihandler < handle

properties
     
    javapath = '';
    authservice = [];
    accessToken = [];
    lastresponse='';
    secretfile = 'AuthKeys.enc';
    apiurl = 'https://saccade.neurosurgery.uiowa.edu/';
    isconnected = false;
    
end
properties (Hidden = true, Access = private )
   % Consumer_Key='';
   % Consumer_Secret='';
   % OAuthToken='';
   % OAuthToken_Secret='';
   keydata = []; % password protected key information
   keyheader = [];
    cryptconfirm = int8(1:8);
    salt = 2015;
end
methods
    
    function me = mwapihandler
        
        
       
        me.javapath = fullfile(pwd,'javaclass','scribe-1.3.7.jar');
        
       import org.scribe.builder.*
       import org.scribe.model.*
       import org.scribe.oauth.*            
      
    end
    
    
        
    
    
    function create_secret(me,keys)
  
        keylabels = {'Consumer_Key','Consumer_Secret'};%,'OAuthToken','OAuthToken_Secret'};
        if nargin < 1 
            keys = inputdlg(keylabels,'modal',false)';
        end
        filename = 'AuthKeys.enc';

        data = cat(1,keylabels,keys);
  
         plaintext = sprintf('\n%s = %s',data{:});
    
      %  plaintext = [me.cryptconfirm(:);int8(double(str(:))-128)];
        encrypt(plaintext,filename);

        me.secretfile = filename;
        
    end
    
    
    function connect(me)%,keys)
        
        %%% Connect to Mediawiki API, loading key from encrypted file or
        %%% from the 2nd argument, a struct with fields that are the key
        %%% labels and values that are the unencrypted keys.
        javaaddpath(me.javapath);

    
        import org.scribe.builder.*
        import org.scribe.model.*
        import org.scribe.oauth.*                

%          keylabels = {'Consumer_Key',...
%                       'Consumer_Secret',...
%                       'OAuthToken',...
%                       'OAuthToken_Secret'};
%            
       
%             if nargin < 2

        if  isempty(me.keydata)
            % header = [];
            [plaintext,me.keydata,me.keyheader,me.secretfile] = decrypt(me.secretfile);
        else

             plaintext = decrypt(me.keydata,me.keyheader);
        end

        fcont = regexp(plaintext,'([^\s]*)\s*=\s*([^\s]*)','tokens');
        fcont = cat(1,fcont{:})';
        keys = struct(fcont{:});
%             end
            
       [res,err] = urlread([me.LabwikiApi,'Special:OAuth/authorize']

        service = org.scribe.builder.ServiceBuilder;
        service = service.provider(api.LabwikiApi);
        service = service.apiKey(keys.Consumer_Key);
        service = service.apiSecret(keys.Consumer_Secret);
        me.authservice =  service.build();
%        me.accessToken = Token(keys.OAuthToken,keys.OAuthToken_Secret);
       
     end
     function respout = get(me,getwhat)
         
         % Send a basic get request
         
          url = [me.apiurl,getwhat];
  
         req = me.request(url,'get');
         
         response = req.send();
            me.lastresponse = response;
            respout = char(response.getBody);
        

     end
     function respout = post(me,getwhat,postdata)
         
         % Send post data.
         % postdata is a struct with optional fields .queryParams, .bodyParams
         % .payload, .header and .payload
         
           url = [me.apiurl,getwhat];
  
          req = me.request(url,'post');
         if nargin > 2
             flds = fieldnames(postdata);
             for k = 1:length(flds)
                  if isempty(postdata.(flds{k}))
                       continue
                  end            
                 
                   switch lower(flds{k})
                       case 'queryparams'
                           javacall = @(a,b)req.addQuerystringParameter(a,b);
                           subfld = fieldnames(postdata.(flds{k}));
                      case 'bodyparams'
                           javacall = @(a,b)req.addBodyParameter(a,b);
                           subfld = fieldnames(postdata.(flds{k}));
                       case {'payload','body'}
                           req.addPayload(postdata.(flds{k}));
                           subfld = {};
                       case 'header'
                            javacall = @(a,b)req.addHeader(a,b);
                           subfld = fieldnames(postdata.(fld{k}));
                       otherwise
                           error('Unrecognized postdata field')                 
                   end
                   
                   for kk = 1:length(subfld)
                      javacall(subfld{kk},postdata.(flds{k}).(subfld{kk})) 
                   end
                   
             end
         end
         %%
         resp = req.send();
         me.lastresponse = resp;
         respout = xml2struct(char(resp.getBody));
              
     end
     
    function req = request(me,url,type)

        % Send a get or post request
        
        import org.scribe.builder.*
        import org.scribe.model.*
        import org.scribe.oauth.*                   
        
        if nargin < 2
            type = 'get';
%         else
%             type = 'post';
        end
        
        if isempty(me.accessToken)
            fprintf('\nNeed to connect to the API first...')
            me.connect;
        end
%             url = [me.apiurl,getwhat];
            switch type
                case 'get'
                    req = OAuthRequest(Verb.GET,java.lang.String(url));
                case 'post'
                    req = OAuthRequest(Verb.POST,java.lang.String(url));
%                     for k = 1:size(params,1)
%                         request.addQuerystringParameter(params{k,1},params{k,2});
%                     end
            end

            me.authservice.signRequest(me.accessToken,req);
%             response = request.send();
%             me.lastresponse = response;
%             respout = char(response.getBody);
                 
        
        
    end


end
end
