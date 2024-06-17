
classdef kstore < handle

%    
%       ks = kstore 
%
% Creates a kstore serial data storage object which allows saving matlab variables 
% to a file in the format described below (which is not compatible with the
% mat-file format). It was developed to allow more rapid serial reading and 
% writing of matlab variables without the unneccessary overhead caused by
% repeated calls to 'save --append' when serially adding or loading multiple variables
% to existing mat files. The  code needed to read and/or write the file is
% prepended to the beginning of every data file, so if you have the data, you
% have the means to access it. Simply copy the code between the hash delimited
% text in the data file to a file called kstore.m if you don't already have
% the file or have an incompatible version. 
%
% How to use this:
%
% ks = kstore('filename.kov'), opens an existing file named filename.kov or
% creates a new one. 
%
% To append a single variable to the file, use the 'append' method:
%
%        ks.append(my_variable1)
%
%        ks.append(my_variable1, newname)
%
% saves the variable with a different name. 
% 
% and 'save' for multiple variables:
%
%       ks.save('var1','var2',...)
% while
%       ks.save    
% 
%  alone saves all variables in the workspace.  
% 
% Note that unlike 'save', 'append(.)' takes the variable
% itself as the argument rather than a string containing its name.
% 
% For the sake of simplicity and efficiency, data can only be appended to
% the file. Saved variables are never overwritten. To replace an old file, 
% delete it first. 
%
% To load one or more variables, use the 'load' method:
%
%        ks.load('my_variable')
% while
%        ks.load 
% 
%  alone loads all variables in the file. 
% 
% The following methods are described for completenes.
% 
% To read data serially from the file into a new variable, use the 'read' method. 
%
%         data = ks.read();
%
% To skip to the next variable without reading it, use 'advance', 
%
%         ks.advance();
%  
% Information about the variables in the file is contained in the struct array
% ks.variables, and for the variable at the current file position in  
% ks.curr_var_data. 
%
% This structure includes the following fields:
% 
%     .name: variable name
%     .class: variable class
%     .size:  size of the variable. Size for struct arrays includes
%             an extra dimension containing the number of fields.
%     .position: file position
%     .flags:   integer whos 1st 3 bits indicate: bit 1 high for nested variables
%               (that is, cell array elements or struct fields), bit 2 for sparse,
%               bit 3 for complex.
%           NOTE: kstore does not currently work with sparse or complex variables.                     
%     .nested: 1 if the current variable is a struct field or a
%            cell array element (also indicated by the flag). 
%
% 
%
%  Additional methods:
%
%    .close():  closes the file. 
%    .scan(): scans variable information in the file. 
%    .start_of_data(): return to the beginning of the data.
%
%   Properties:
%
%    .filename: file name
%    .path: file path
%    .isopen: true when a file is open for reading and/or writing
%    .access: current read/write access: 'r+' and 'w+'- read/write, 'r' -
%             read only
%    .file_position: current file position
%    .eof: true if file position is at the end of file. 
%    .current_variable: name of the variable or fieldname for struct data
%                at the current position.
%    .variable_index: index of the variable in the file. 
%    .read_precision: numeric data are read into matlab with this precision.
%                 Options are 'double','single' or any int class and 'same'
%               to readvariables in the same format they are written to file. 
%               
%  Kstore will only save variables of the standard data classes (double,
%  single, intN, char, etc) as well as cell and struct. It will NOT save
%  function handles or other non-standard classes, and will skip over them.
%
%
%  The file format is explained below. 
% 
%
%%%%% File Format %%%%%%
%
% File Header: 
%       
%       Byte offset to the end of the file header       uint32 x 3 
%             thrice repeated    
%       Contents of this m-file                         uint8  x Number of
%                                                               bytes
% File Body: 
%   
%     Variable Header:                                 
% 
%       Start code:     'zscmk'                         uchar x 5
%       Name length: length of name thrice              uint8 x 3
%       Name:  variabe name                             uchar x length
%       Byte offset to the end of the variable header   uint32 x 3 
%             thrice repeated                       
%       Type:   variable class                          uchar x 7
%       Flags:   flags for variable attributes          uchar x 1
%              bits 1-3 indicate 1-nested, 2-sparse, 3-complex
%       Ndim:   number of dimensions                    uint64 x 3
%               repeated thrice.
%       Size:   number of elements in each dimension    uint64 x Ndim 
%       Padding:  unused header space (or other data)   variable 
%       
%       End code:       'kmcsz'                         uchar x 5
%       Offset: number of bytes in the variable body    uint64 x 3
%               thrice repeated.
%     Variable Body:                                    Bytes specified in 
%                                                       variable header
  

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/misctools/kstore.m $
% $Revision: 696 $
% $Date: 2016-01-21 21:08:55 -0600 (Thu, 21 Jan 2016) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2013

    properties ( Dependent=true )
        filename = '';
        path = '';
       isopen = [];
       access = [];
       file_position = [];
       eof = [];
       read_precision = [];
    end
    
    properties (Dependent = true, Hidden = true)
       header=[]; 
       fileheader = [];       
    end
    properties ( Hidden = true, SetAccess = private )
        accessval =0;  % 0-closed, 1-read, 2-write,3-read/write
        fid = -1;
        newfile = false;
        mcode_preamb = '###### COPY THIS MATLAB CODE TO READ/WRITE THIS FILE -- DO NOT EDIT ANYTHING HERE! ######';
        data_preamb = '###### BEGIN DATA ######';
        start_code='zscmk'; % This is included at the start of every variable header
                              % and the reverse at the end of every header.
       precval = 'double';
       scanned = false;
       prev_file = [];
       
       uuid = [];
       writeblockfile = '';
    end
    
    properties ( SetAccess = private )
       
        current_variable = '';
        variables = struct('name',[],'name_length',[],'class',[],'flags',[],...
                            'dimensions',[],'size',[],'position',[],'nested',false);
        variable_index=0;
        curr_var_data = struct('name',[],'name_length',[],'class',[],'flags',[],...
                            'dimensions',[],'size',[],'position',[],'nested',false);
     
    end
    
    methods
        function me = kstore(varargin)
            
            if usejava('jvm')
                me.uuid = char(java.util.UUID.randomUUID);
            else
                me.uuid =  num2str(matlab.internal.timing.timing('cpucount'));
            end

            if nargin < 1 || isempty(varargin{1})
               fn  = uifileselect({'kov','Kovach Store Files (*.kov)'}); 
            else
                fn = varargin{1};
            end
            
            if nargin < 2 || isempty(varargin{2})
                perm = 'r';
            else
                perm = varargin{2};
            end
            
            [pth,fn,ext]=fileparts(fn);
            if isempty(ext)
                ext = '.kov';
            end
            if isempty(pth)
                pth = cd;
            end
            filename =[pth,filesep,fn,ext];
            me.newfile = ~exist(filename,'file');
            me.writeblockfile = fullfile(pth,sprintf('.%s%s',fn,ext));
            if me.newfile
                me.openfile(filename,'w+');
                me.variables = me.variables([]);
                me.curr_var_data= me.variables([]);
            else
                me.openfile(filename,perm);
                me.scan;
            end            
      
        end
        %%%%%
        
        function openfile(me,fn,varargin)
            
            if me.isopen > 0
                fprintf('File already open');
                return
            end
            
            if nargin < 2
                fn = me.prev_file;
            end
            if nargin < 3
                access = 'r';
            else
                access = varargin{1};
            end
                    
            
            me.fid = fopen(fn,access);
            
            if me.fid < 0             
                error('Failed to open file.')
            end
            
            switch me.access
                case {'w','w+','wb+'}
                    % First get mfile contents to save with the file
                    % 'cause why not?
                    
                    if isequal(me.access,'w') || me.newfile  
                        mfn = which(mfilename);
                        mfid = fopen(mfn,'r');
                        mfcont = fread(mfid,'uchar=>uchar');
                        fclose(mfid);

                        nb = numel(mfcont);
                        me.writeblock;
                        fwrite(me.fid,nb*[1 1 1],'uint32');
                        fprintf(me.fid,'\n%s\n',me.mcode_preamb);
                        fwrite(me.fid,mfcont,'uchar');
                        fprintf(me.fid,'\n%s\n',me.data_preamb);
                        fwrite(me.fid,[0 0 0],'uint64');
                        me.unblock;
                    end
            end
            me.prev_file = fn;
            
        end
        
        
        
        %%%%%%%
        function save(me,varargin)
            %%% This behaves like the matlab save function. Input is
            %%% a variable name as a string. With no input the caller's
            %%% workspace is retrieved.
            
            if nargin < 2
                vars = evalin('caller','who');
            else
                vars = varargin;
            end
            for i = 1:length(vars)                
                me.append(evalin('caller',vars{i}),vars{i});
            end
        end
        
        %%%%%%%
        function append(me,data,varname,nested)
           
%             data = stripfunctions(data);
            if isequal(me.access ,'r') || isequal(me.access,'rb') 
               
                inp = questdlg('Are you sure you want to modify this file (set .access to ''r+'' to avoid this question)?');
                if strcmp(inp,'Yes')
                    me.access = 'r+';
                end
                
            end
            
            if isnumeric(data)
                type = 'numeric';
            else
                type = class(data);
            end
            if nargin < 3
               varname = inputname(2);
            end
           
            if nargin < 4  %%% Nested is true of structure fields and cells.
                           %%% These should not be included in the variable
                           %%% list. 
                nested = 0;
            end
            
          
            
            mdat = whos('data');
            
            flags = nested + 2*mdat.sparse + 4*mdat.complex;
            if flags > 1
                warning('Not yet able to handle sparse or complex data...coming soon.')
                fprintf('\nSkipping...')
                data = [];
                mdat = whos('data');
                flags = nested + 2*mdat.sparse + 4*mdat.complex;
            end
            
            cl = class(data); 
            cl(8:end)='';
            cl(end+1:7)='_';
            

            vardata = me.curr_var_data([]);
            
            
            vardata(1).name = varname;
            vardata.name_length = numel(varname);
            vardata.class = cl;
            vardata.flags = flags;
            vardata.dimensions = ndims(data) + isstruct(data); %for structs, additional dimension for field number
            vardata.size = size(data);
            vardata.nested = nested;
            if isstruct(data)
               vardata.size(end+1) = length(fieldnames(data));
            end
            
            hd = me.header;
            
%             hd = me.header;
            
            if ~me.eof
                fseek(me.fid,0,1);
            end
            
            unblock = true; 
            me.writeblock;
            
            vardata.position = ftell(me.fid);
            fwrite(me.fid,me.start_code,'uint8');
            
                        
            %%% Write the header using the template
            for i = 2:length(hd)-2                
                to = hd{i,3};
                xout = vardata.(hd{i,1});
                fwrite(me.fid,xout*ones(size(hd{i,2})), to );                
            end
            fwrite(me.fid,me.start_code(end:-1:1),hd{end-1,3});
            
            
            
            flp0 = me.file_position;
            %%%%%%% place holder for number of bytes                      
            fwrite(me.fid,zeros(size(me.header{end,2})),'uint64');  
        
            %%%%%%% now write the data
             switch type
                
               case 'logical'

                    fclass = 'ubit1';
                    fwrite(me.fid,data,fclass);
                    
                case 'char'

                    fclass = 'uchar';
                    fwrite(me.fid,data,fclass);
                    
%                  case 'function_handle'
%                      fclass = 'uchar';
%                      me.append(func2str(data),'funcstr',1);
                     
                case 'numeric'

                    fclass = mdat.class;
                    fwrite(me.fid,data,fclass);

                case {'cell','struct'} % For cell and struct arrays each field or element is treated like a separate variable
                      
                      if isstruct(data)
                          
                           fldn = fieldnames(data);
                          for k = 1:numel(data)
                                for i = 1:length(fldn)
                                     me.append(data(k).(fldn{i}),sprintf('%s',fldn{i}),1);
                                 end
                           end
                       elseif iscell(data)
                           for i = 1:numel(data)
                                me.append(data{i},sprintf('%s',varname),1);                      
                           end
                      end
                   unblock = false;     %#ok<*NASGU>
%                  case 'function_handle'
%                   warning('\nSkipping variables of class %s. The variable ''%s'' will NOT be stored.',type,varname) %#ok<WNTAG>
%                    fclass = 'uchar';
%                     fwrite(me.fid,' ',fclass);
                     
                 otherwise 
         
                    beep
                    warning('\nSkipping variables of class %s. The variable ''%s'' will NOT be stored.',type,varname) %#ok<WNTAG>
%                     fwrite(me.fid,' ','uchar');
%                     sz = [0 0];
%                     fseek(me.fid,vardata.position,-1);
%                     return
             end
             
             %%%%% Write the variable size in bytes at the start and return
             %%%%% to the current position
            flp = me.file_position;
            hd = me.header;
            nb = flp -flp0 - sum(hd{end,2})*me.nbytes(hd{end,3});
            fseek(me.fid,flp0,-1);
            fwrite(me.fid,nb*me.header{end,2},'uint64');
            fseek(me.fid,flp,-1);
            
            %%%% Update variable information
            me.current_variable=varname;
            me.curr_var_data=vardata;
            if ~nested
                me.variable_index = me.variable_index+1;
                me.variables(me.variable_index)=vardata;
            end  
            
            if unblock
                me.unblock;       
            end
        end
     
        %%%%%%%%%
        
        function close(me)
            % Close the file
            if ~isempty(me.access)
                fclose(me.fid);
            end                   
        end
        
        %%%%%%%%%
        function  ex = get.isopen(me)
            % Close the file
            fn = fopen(me.fid);
            
            ex= ~isempty(fn);
                
        end
        %%%%%%%%%%%%%
        function  perm = get.access(me)
            % Retrieve access permission
           
              [~,perm] = fopen(me.fid);
           
        end
        %%%%%%%
        function a = get.eof(me)
            if me.isopen
                a = feof(me.fid);
            else
                a = [];
            end
        end
        
        %%%
        function   set.access(me,acc)
            % Reset access permission
            [fn,perm] = fopen(me.fid);
            
            if strcmp(acc,'w')
                warning('Opening as read/write instead of write only.') %#ok<WNTAG>
                acc = 'r+';
            end
            
            if ~strcmp(acc,perm)
                fclose(me.fid);
                try
                    me.fid = fopen(fn,acc);
                    if me.fid < -1
                        error('Failed to open')
                    end
                    me.start_of_data();
                catch err
                    me.fid = fopen(fn,perm);
                    rethrow(err)
                end
            end
                                      
        end
        
        function  file = get.filename(me)
            fn = fopen(me.fid);
            
            [~,file,ext] = fileparts(fn);
            
            file = [file,ext];
            
        end
        function  pth = get.path(me)
            fn = fopen(me.fid);
            
            [pth,~,~] = fileparts(fn);
            
        end
            
        function a = get.read_precision(me)

            a= me.precval;
                  
        end
        function set.read_precision(me,a)
            
            try
                
                realmax(a); %to check that this is a valid precision on this system
                  
            catch err1 %#ok<NASGU>
                try 
                    intmax(a)
                catch err2
                    error('%s is not a recognized precision (e.g. try ''double'' or ''single'')',a)
                end
             end
                me.precval = a; %Read all numeric value into matlab with this precision.
        end
        
        function  pos = get.file_position(me)

            if ~isempty(me.access)
               pos = ftell(me.fid); 
            else
                pos = [];
            end
            
        end
        
        %%%%%
        function delete(me)
            
            me.close();
            
        end
     
        
        %%%%%
            
        function  start_of_data(me)
           
            frewind(me.fid)
            nb = fread(me.fid,3,'uint32');
            
            unb = me.vote(nb);

            [~] = fgets(me.fid);
            [~] = fgets(me.fid);
            stat =fseek(me.fid,unb,0);
            if stat < 0 
                ferror(me.fid);
            else
                [~] = fgets(me.fid);
                [~] = fgets(me.fid);
            end
            me.variable_index = 0;
            me.current_variable = '';
            me.curr_var_data = me.curr_var_data([]);
            me.advance;
            
        end

        
        %%%%%%%
        function  advance(me,offset)
            
            
            %%% Advance by offset and read the next variable header. 
            %%%
            
            if me.file_position == 0
                me.start_of_data();                
            end
            if nargin < 2
                if ~me.eof  %%% with no argument, assume that the offset is 
                            %%% written at the current file location.
                    offset = me.vote(fread(me.fid,3,me.header{end,3}));
                else
                    offset = 0;
                end
            end           
            
            fseek(me.fid,offset,0);
            
            vardata = me.variables([]);
            %%%%% Advance to the next variable 
            vardata(1).position = me.file_position;
            stc = fread(me.fid,5,'uchar=>char')';
             if isempty(stc)
%                 fprintf('\nAt end of file. No further data to read.\n')
                return
             elseif ~strcmp(stc,me.start_code)
                error('Error reading from file... either I got distracted and lost my place or the file is corrupt.');
            end
            
            hd =  me.header;
            xin = [];
        
            for i = 2:length(hd)-2
                if hd{i,2} == 0
                    n = me.vote(double(xin));
                else
                    n = sum(hd{i,2});
                end
                
                from = hd{i,3};
                if strcmp(from,'uchar')                    
                    to = 'char';
                else
                    to = me.read_precision;
                end
                xin = fread(me.fid, n,[from,'=>',to])';
                vardata.(hd{i,1}) = xin;
            end
            
            vardata.nested = bitand(vardata.flags,1);
            vardata.class = vardata.class(vardata.class~='_');
            
            endcode = fread(me.fid, length(me.start_code),'uchar=>char')';        %#ok<NASGU>
            
       
            me.curr_var_data = vardata;
            me.current_variable = vardata.name;
            if ~vardata.nested
                me.variable_index = me.variable_index+1;
            end
        end
        
        %%%%
        function scan(me)
            
            me.start_of_data;
            while ~me.eof
                if ~me.eof && ~me.curr_var_data.nested
                    me.variables(me.variable_index) = me.curr_var_data;
                end
                me.advance;
            end
            me.scanned = true;
            me.start_of_data;
            
        end
        %%%%%%%%%
        function rewind(me)
           
            frewind(me.fid);
        end
        
        %%%%%%%%%%  Read the current variable and return it
        
        function data = read(me)
            
            if me.eof
                fprintf('\nat the end of the file...')
                data = [];
                return
            elseif me.variable_index == 0
                me.advance;
            end
            vardata = me.curr_var_data;
            fclass = vardata.class(vardata.class~='_');
            
            sz = vardata.size;
            n = prod(sz); 
            hd = me.header;
            fseek(me.fid,sum(me.header{end,2})*me.nbytes(hd{end,3}),0); %%% Advance skipping to the beginning of the data
            
            switch fclass
            
                case 'struct'
                   
                   fldn = sz(end);
                   sz = sz(1:end-1);
                   n = prod(sz); 
                   for k = 1:n     
                      for i = 1:fldn
                          me.advance(0);     
                           svardata = me.curr_var_data;
                           data(k).(svardata.name) = me.read; %#ok<AGROW>
                       end 
                   end
                   if ~exist('data','var')
                       data = [];
                   end
%                 case 'functio'
%                     svardata = me.curr_var_data;
%                     me.advance(0);
%                     data = me.read;
%                     data = str2func(data);
%                     me.advance(0);
                case 'cell'
                    
                   
                   data = cell(sz);
                   for k = 1:n                   
                       me.advance(0);         
                       data{k} = me.read;
                       
                   end
                case 'char'
                      data = fread(me.fid,n,'uchar=>char');
                case 'logical'
                      data = fread(me.fid,n,'ubit1=>logical');
%                 case 'functio' % Currently functions are not stored
%                       data = fread(me.fid,1,'uchar=>char');
%                       sz = size(data);
                otherwise
                    try
                      data = fread(me.fid,n,[fclass,'=>',me.precval]);
                    catch err %#ok<NASGU>
                        data = [];
                        sz = [0 0];
%                         sz = [1 1];
                    end
            end
            
            if ~isa(data,'function_handle')
            data = reshape(data,double(sz));
            end
            if ~vardata.nested
                me.advance(0);     
            end            
            
        end
          
        
        %%%%%%%% Load data in the workspace like the matlab load function
        function varargout = load(me,varargin)
           
            if ~me.scanned
                me.scan;
            end
            
            if nargin < 2
                loadvars = {me.variables.name};                
            else
                loadvars = varargin;
            end
            
            varnames = {me.variables.name};
            getvar = ismember(varnames,loadvars);
            
            if sum(getvar) < length(loadvars)
                missing = setdiff(loadvars,varnames(getvar));
                warning('The following variables were not found in this file: %s',...
                        sprintf('\n%s',missing{:})); %#ok<WNTAG>
            end
            
            %%% emulate behavior of load function 
            for i = find(getvar)
                var = me.variables(i);
                fseek(me.fid,var.position,-1);
                me.advance(0);
                me.variable_index = i;
                data = me.read;
                if nargout > 0
                   varargout{1}.(var.name) = data;
                else
                    assignin('caller',var.name,data);
                end
                
            end
            
        end
        
          %%%%%%%% Go to start of specified variable
        function  gotovar(me,var)
           
            if ~me.scanned
                me.scan;
            end
            
          
            
            varnames = {me.variables.name};
            getvar = ismember(varnames,var);
            
            if sum(getvar) ==0
                warning('The following variables were not found in this file: %s',...
                        sprintf('\n%s',var)); %#ok<WNTAG>
            end
            
            %%% advance to the variable 
            var = me.variables(getvar);
            fseek(me.fid,var.position,-1);
            me.advance(0);
            me.variable_index = find(getvar);

       
            
        end  
        
        
        %%%%%%
        function hd = get.header(me)
            
            %%%% Defines the variable header
            hd = {'start_code',length(me.start_code),'uchar'
                'name_length',[1 1 1],'uint8'
                'name',0, 'uchar'
                'class',7,'uchar'
                'flags',1,'uint8'
                'dimensions',[1 1 1],'uint64'
                'size',0,'uint64'
                'end_code',length(me.start_code),'uchar'
                'bytes_in_data',[1 1 1],'uint64'};
        end
        
         %%%%%%
        function fhd = get.fileheader(me) %#ok<MANU>
            
            
            %%%% Defines the file header
            fhd = {'byte_offset',[1 1 1],'uint32'
                   'text',0,'uchar'};
                
        end
        
        %%%%% 
        function writeblock(me)
            % Deny any asynchronous kstore process write access.
            while ~me.checkblock
            end
        
            ffid = fopen(me.writeblockfile,'w');
            fprintf(ffid,me.uuid)';
            fclose(ffid);
            
        end
        %%%%%%%
        function unblock(me)
           % Allow write access to other processes
           while ~me.checkblock
           end
           
           delete(me.writeblockfile)
          
        end
        %%%%%%
        function isfree = checkblock(me)

            agelimit = 600; % Assume a process that has been blocking for longer than this many seconds is crashed.

            
            dd = dir(me.writeblockfile);
            isfree = isempty(dd);
            if ~isfree
                ffid = fopen(me.writeblockfile,'r');
                proc = fread(ffid,'uchar=>char')';
                fclose(ffid);
                
                isfree = strcmp(proc,me.uuid) || ((now-dd.datenum)*3600*24 < agelimit);
            end
                
        
        end
            %%%%%%%
    end

    
   %%%%%%%% static methods
   methods ( Hidden = true, Access = private, Static = true )
 
    
        function unb = vote(nb)
            % C
            unb = unique(nb);
            if length(unb) == 2
                warning('A corrupted value for number of header bytes has been detected.') %#ok<WNTAG>
                unb = mode(nb);
            elseif length(unb) == 3
                error('Number of bytes in header is corrupted. Halting.')
            end
        end
        
        function nb = nbytes(class)
            
            switch class
                case {'uchar','char','uint8','int8'}
                    nb = 1;
                case {'uint16','int16'}
                    nb = 2;
                case {'uint32','int32','single'}
                    nb = 4;
                case {'uint64','int64','double'}
                    nb = 8;
                otherwise
                    error('unknown type')
            end
        end
                    
   end
   
end
            