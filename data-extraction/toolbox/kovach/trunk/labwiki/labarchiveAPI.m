
classdef labarchiveAPI < handle

    % Interface for the labarchive API. 
    %
    % USAGE:
    %
    % Initialilze a labarchiveAPI object with:
    %   lapi = labarchiveAPI;
    % 
    % Enter the decryption password (try the usual password for Lab)
    % 
    % Then give your labarchive login credentials.  
    % One cautionary note: the labarchive API transfers credential information in a way 
    % that isn’t very secure, using a GET call, meaning that the information is passed as
    % within the requested url.  The reason this is bad practice is that servers often 
    % write the url to a log file, meaning that after this step their servers 
    % may contain an unencrypted text file with your credential information.  So if you 
    % do use this, I would recommend changing your labarchive password to one that you don’t
    % use for a lot of other things, and especially don’t use a password that is the same or 
    % similar to the one for the email address you use to log in to labarchive as, otherwise, anyone 
    % who breaks into their servers might be able to hack into your email account right away.
    % To avoid repeatedly sending credentials in this insecure way, an encrypted file is written the first 
    % time you log in that contains the unique encryption key needed to
    % connect to the server.
    %
    % Once you are logged in you can choose an action from the list under “API Classes” in the link above and specify the parameters ,except ‘uid' and ‘sig', which the interface handles for you. For example to retrieve the basic  information for 405-002, the call would be
    % Specify the action 
    %    lapi.current_method=‘entries/entry_info’;
    %
    % Code for the entry you want to retrieve
    %    entry_identifier_code='NTQ0MS44fDIyODE4NC80MTg2L0VudHJ5UGFydC8yNDU2MTc2N3wxMzgxMy44’;
    %
    % Create the query
    %    lapi.query_params =struct('eid’,entry_identifier_code,'entry_data','true');
    %
    % Generate the url for the api call
    %    lapi.create_url; 
    %
    % Append an encrypted signature to the url
    %    lapi.sign_url;   
    %
    % Get the result
    %    xml_result= webread(lapi.url);  % result from the labarchive server
    %    parsed_result = lapi.parse_entry(xml_result); % Parse the result and return a more useful struct output
    % 
    % The information for the form identifier code above is then contained as fields in
    %    parsed_result.parsed
    % 
    % The entry  identifier codes you can get by populating a tree containing all entries in the following way:
    % 
    % Populate the list of notebook entries with
    %   lapi.get_notebook(1,0);
    %
    % Pick the subject you want based on lapi.notebooks(1).tree.label
    % Get the full contents of the kth tree (containing the subject you want)
    %   tree = lapi.get_tree(1,Inf,k);
    %
    % Pick a block based on tree.label — say, entry 10
    %   entry_identifier=tree(10).entry(1).entryID
    %
    
    %C. Kovach 2018
    properties
        baseurl = 'api.labarchives.com';  %API access site
        akid='uiowa_kovach';         %API identifier
        keyfile = 'labarchive.enc';  %This contains the api secret key, encrypted
        userfile = 'lapi_users.enc'; %After the initial login, user-identifying information will 
                                     %be stored encryped in this file so full credentials don't 
                                     %have to be resent at every login
        userfilepath='';       %Path to above file. Defaults to that of labarchiveAPI.m     
        user = '';             %User name
        userid = '';           
        current_method = '';   %Currently specified action. These are described under "API Classes" here: https://mynotebook.labarchives.com/share/LabArchives%20API/NS4yfDI3LzQvVHJlZU5vZGUvMTF8MTMuMg==
        query_params;          %Parameters required by the above action, specified as a struct with query_params.(parameter_name)=parameter_value. Alternatively can be a cell with {'parameter_name','parameter_value',...}
        url='';                %Current url
        signed=false;          %Whether the current url is signed
        notebooks = [];        %Notebooks associated with the active account
        logged_in = false;     %Whether successfully logged in
        tree=[];               
        errmsg='';             %Last caught error 
        active_notebook=1;     %Index of the active notebook
    end
    properties (Dependent = true)
       timestamp 
       
    end
    properties (Hidden = true)
        
        mac = '';
    
    end
        
    properties (Access = private)
       salt = 'tp39dc48b3_a0a1_4cd7_8271_5b181ed476d5';
    end
    
    methods  
        
        function me = labarchiveAPI
            try
              me.createMAC();
              if isempty(me.userfilepath)
                 me.userfilepath = fileparts(which(mfilename));          
              end
              
              me.login();
              me.get_notebooks(1:length(me.notebooks),0);
            catch err
                me.errmsg = err;
                warning(err.message)
            end
        end
        function createMAC(me)
            import java.security.InvalidKeyException;
            import java.security.NoSuchAlgorithmException;
            import java.security.SignatureException;
            import java.util.Formatter;

            import javax.crypto.Mac;
            import javax.crypto.spec.SecretKeySpec;
            
            
            
            sec = me.getsec;
            seckey = SecretKeySpec(uint8(sec), 'HmacSHA1');
            me.mac = Mac.getInstance('HmacSHA1');
            me.mac.init(seckey);
            
            
        end
        function get_notebooks(me,notebook,depth)
            if nargin < 3 || isempty(depth)
                depth = 0;
            end
            if nargin < 2 || isempty(notebook)
                notebook = 1:length(me.notebooks);
            end
            if ~isnumeric(notebook)
                notebook = find(ismember({me.notebooks.id},{notebook.id}));
            end
            for k = 1:length(notebook)
                nbtree = me.get_tree(notebook(k),depth);               
                me.notebooks(notebook(k)).tree = nbtree;
            end
        end
        function login(me,recredential)
             
            if nargin < 2
                recredential = false;
            end
            
            if isempty(me.mac)
                me.createMAC;
            end
            persistent uid nbks loginid
            
            xml='';
           if isempty(uid) && isempty(me.userid) || recredential
               
             
           
                newloginid = input(sprintf('Enter LabArchive login [%s]:',loginid),'s');
                if ~isempty(newloginid)
                    loginid = newloginid;
                end
                pw = deblank(getpw('Enter LabArchive Password:'));
                if exist(fullfile(me.userfilepath,me.userfile),'file') && ~recredential % Store the uid to an encrypted file to avoid re-sending login credentials in insecure GET  
                    userdat = me.getsec(fullfile(me.userfilepath,me.userfile),me.akid);
                    re = regexp(userdat,sprintf('%s:(?<password>.*?)::::(?<header>[^\n]*)',loginid),'names','once');
                else
                    re=[];
                    userdat = '';
                end
                if isempty(re) 
                    me.current_method='users/user_access_info';
                    me.query_params.login_or_email = loginid;
                    me.create_url;
             	    me.sign_url;
                    me.url = sprintf('%s&password=%s',me.url,pw);
                   
                    xml = webread(me.url);
                    uid = regexp(xml,'<id>([^<]*)</id>','tokens','once');
                    uid = uid{1};
                    
                    if ~isempty(me.userfile)
                        [uicphr,uihdr] = encrypt(uid,'',[me.salt,pw]);
                        if ~isempty(regexp(userdat,sprintf('%s:',loginid),'once'))
                            userdat = regexprep(userdat,sprintf('%s:[^\n]*',loginid), sprintf('%s:%s::::%s',loginid,me.base64enc(uicphr),me.base64enc(uihdr)));
                        else
                            userdat = sprintf('%s\n%s:%s::::%s',userdat,loginid,me.base64enc(uicphr),me.base64enc(uihdr));
                        end
                        me.create_secret_file(fullfile(me.userfilepath,me.userfile),userdat,me.akid);
                    end
                else
                    uid = decrypt(me.base64dec(re.password)',me.base64dec(re.header)',[me.salt,pw]);
                    me.userid = uid;
                    
                end
            elseif isempty(me.uid)
                me.userid = uid;            
           end
           if isempty(xml)
                me.current_method='users/user_info_via_id';
                me.query_params=struct('uid',uid);
                me.create_url;
                me.sign_url;
                xml = webread(me.url);
           end
            me.userid=uid;
            fullname=   regexp(xml,'<fullname>(.*)</fullname>','tokens','once');
            nbks=   regexp(xml,'<notebooks(.*)</notebooks>','tokens','once');
            nbks = regexp(nbks{1},'<notebook>(.*?)</notebook>','tokens');
            nbks = regexp([nbks{:}],'<id>(?<id>.*?)</id>.*<name>(?<name>.*?)</name>','names');
            nbks=  [nbks{:}];
            me.notebooks = nbks;
            me.logged_in=true;
            me.user= fullname{1};
            
        end
        
        function create_url(me)
            
          me.url = sprintf('https://%s/api/%s?%s',me.baseurl,me.current_method,me.create_query);
          me.signed=false;
        end
        
        function sign_url(me)
           
            expires = me.timestamp;
            [~,method] = strtok(me.current_method,'/');
            
            if ~isempty(me.userid)
                uidstr = sprintf('uid=%s&',me.userid);
            else
                uidstr='';
            end
            sig = me.base64enc(me.mac.doFinal(uint8([me.akid,regexprep(method,'/',''),expires])));
            me.url = sprintf('%s&akid=%s&expires=%s&%ssig=%s',me.url,me.akid,expires,uidstr,sig);
            me.signed = true;
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
                    qstr=sprintf('%s%s=%s&',qstr,fldn{k},me.query_params.(fldn{k}));
                end
                qstr(end)='';
            elseif iscell(me.query_params)

                qstr = sprintf('%s=%s&',me.query_params{:});
                qstr(end)='';
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
        
        function treeout = get_tree(me,notebook,depth,tree)
            
            if nargin < 2 || isempty(notebook)
                notebook=me.notebooks(me.active_notebook);
            end
            if isnumeric(notebook)
                nbi = notebook;
                notebook=me.notebooks(notebook);
            else
                nbi = find(ismember({me.notebooks.id},{notebook.id}));
            end
            if nargin < 4
                tree.treeID = '0';
                tree.label = 'root';
                tree.isPage='false';
                tree.leafnode = false;
            elseif isnumeric(tree)
                if ~isfield(notebook,'tree')
                    me.get_notebooks(notebook,0);
                end
                tree = me.notebooks(nbi).tree(tree);
            end
            if nargin < 3
                depth = Inf;
            end
            me.current_method = 'tree_tools/get_tree_level'; 
           
         %   for k = 1:length(notebook)
            me.query_params=struct('nbid',notebook.id,'parent_tree_id',tree.treeID);
            me.create_url;
            me.sign_url;
            xml = webread(me.url,weboptions('TimeOut',30));
            re = regexp(xml,'<tree-id.*?>(?<treeID>.*?)</tree-id>.*?<display-text.*?>(?<label>.*?)</display-text>.*?<is-page.*?>(?<isPage>.*?)</is-page>','names');
%            fprintf('\n%s: %s',tree.label,re(kk).label);
            for kk = 1:length(re)
               re(kk).leafnode = strcmp(re(kk).isPage,'true');
               if ~re(kk).leafnode && depth >0
                 re(kk).branch = me.get_tree(notebook,depth-1,re(kk));
               elseif re(kk).leafnode
                   re(kk).entry = me.parse_entry(re(kk),notebook);
               end
            end
            treeout = re; 
           % end
            
        end
        
        function entry = parse_entry(me,tree,notebook)
            
            if ~ischar(tree)
                me.current_method = 'tree_tools/get_entries_for_page'; 

                me.query_params = struct('page_tree_id',tree.treeID,'nbid',notebook.id,'entry_data','true');
                me.create_url;
                me.sign_url;
            else
                xml=tree;
            end
            xml = webread(me.url, weboptions('Timeout',30));
            re = regexp(xml,'<entry>(.*?)</entry>','tokens');
            
            entry= struct('data','','last_modified_by','','modified','','entryID','','entryURL','','parsed','');
            entry=entry([]);
            for k = 1:length(re)
                re2 = regexp(re{k}{1},'<entry-data.*?>(?<data>.*?)</entry-data>.*?<last-modified-by>(?<last_modified_by>.*?)</last-modified-by>.*?<updated-at.*?>(?<modified>.*)</updated-at>.*?<eid>(?<entryID>.*?)</eid>.*?<entry-url>(?<entryURL>.*?)</entry-url>','names');
                if ~isempty(re2)
                    if ~isempty(re2.data)
                        try
                        frmdat = parsejson(re2.data);
                        re2.parsed = frmdat.form_data;
                        catch
                            continue
                        end
                          
                    else
                        re2.parsed=struct([]);
                    end
                    entry(k)=re2;
                end
            end
            
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
