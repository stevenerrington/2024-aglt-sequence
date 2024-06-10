
classdef result < mwapi2
    % This class is for managing results and uploading them to a wiki
    
    % C. Kovach 2018
    
    properties
        
        page_name
        values
        analysis
        file_url
        file_version
        mediafile
        pdf
        pdfpage
        description
        repository_url
        revision_number
        overwrite = false;
        pvalue;
        available_response_types;
        available_subjects;
        available_blocks;
        available_contacts;
    end
    properties(Dependent = true)
        subject
        block
        contacts
        response_type
    end
     properties (SetAccess=protected)
         guih
         handles=struct;
          fields = {'Subject','Block','Contacts','Analysis','Values','Summary PDF','PDF Page','Summary Media','Repository URL','Revision Number','Response Type','Description','Submit','Pvalue'};
        subjid='';        
        blockid = '';
        responsetypeid='';
        contactnums;
     end
     
    methods
        
        function [fig,handles] = gui(me,fields,callback)
           
            
            if ~isempty(me.guih) && ishandle(me.guih) && nargin <2
                figure(me.guih);
                me.fillgui;
                return
            end
            
            fig = figure('units','normalized');
            if nargin < 2
                me.guih = fig;
                callback = @(varargin)me.mwupdate(varargin{:});
                fields = me.fields;
                if isempty(me.available_response_types) 
                    resp = me.askargs({'Category','Response Type'});
                     if ~isempty(resp)
                         me.available_response_types = {resp.main};
                     end
                end
                if isempty(me.available_subjects)
                    subjs= me.askargs({'Category','Subject'},'Has ID',{'sort','Has_ID'});
                    subjs=subjs(end:-1:1);
                    fakesubs =  me.askargs({'Is a fake subject','True'},'Has ID');
                    subjids = setdiff([subjs.Has_ID],[fakesubs.Has_ID],'stable');
                    me.available_subjects = subjids;
                end
            end
            w = 1/(length(fields)+5);
            pos = [.02 1-w*1];
            
            for k = 1:length(fields)
                switch fields{k}
                    case 'Description'
                        dboxheight = 3;
                          uicontrol('style','text','units','normalized','position',[pos+ [0,-.025], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
                         pos = pos + [0 -1]*dboxheight*w;
                         ui = uicontrol('style','edit','units','normalized','position',[pos+ [0,-.025],.9,dboxheight*w],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left','max',2);
                         pos = pos+[0 -dboxheight+1]*w;
                    case 'Summary Media'
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.2,0],.35,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left','Tag',fields{k});
                          uicontrol('style','pushbutton','units','normalized','position',[pos+ [0,0], .2,w*.9],'string',[fields{k},':'],'HorizontalAlignment','left','callback',@(varargin)me.upload('',ui));
                         pos = pos+[0 -1]*w;
                    case 'Summary PDF'
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.2,0],.3,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left','Tag',fields{k});
                          uicontrol('style','pushbutton','units','normalized','position',[pos+ [0,0], .2,w*.9],'string',[fields{k},':'],'HorizontalAlignment','left','callback',@(varargin)me.upload('',ui,'pdf'));
                    case 'PDF Page'
                         uicontrol('style','text','units','normalized','position',[pos+[.45 0], .1,w*.9],'string','Page','HorizontalAlignment','left');
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.5,0],.05,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left');
                         pos = pos+[0 -1]*w;
                    case 'Repository URL'
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.2,0],.3,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left','Tag',fields{k});
                          uicontrol('style','text','units','normalized','position',[pos+ [0,0], .2,w*.9],'string',[fields{k},':'],'HorizontalAlignment','left');
                    case 'Revision Number'
                         uicontrol('style','text','units','normalized','position',[pos+[.45 0], .1,w*.9],'string','Rev.','HorizontalAlignment','left');
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.5,0],.05,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left');
                         pos = pos+[0 -1]*w;
                    case 'Subject'
                        val =find( strcmp(me.subject,me.available_subjects))+1;
                        if isempty(val)
                            val = 1;
                        end
                                 uicontrol('style','text','units','normalized','position',[pos+ [0,-w*.9], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
%                           uicontrol('style','text','units','normalized','position',[pos + [.2,0], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
                         ui = uicontrol('style','popupmenu','units','normalized','position',[pos+[.2,-w*.9],.25,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.subject_callback(a,b),'string',{'Subjects',me.available_subjects{:}},'HorizontalAlignment','left','value',val);
%                          pos = pos+[0 -1]*w;
                     case 'Submit'
                         ui= uicontrol('style','pushbutton','units','normalized','position',[pos+ [.3,0], .3,w*1.5],'string',fields{k},'HorizontalAlignment','left','callback',@(varargin)me.uisubmit);

                     case 'Block'
                         if ~isempty(me.block)
                            val = find(strcmp(me.block,{me.available_blocks.main}));
                         else
                             val=1;
                         end
                          uicontrol('style','text','units','normalized','position',[pos + [.5,0], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
                         if ~isempty(me.available_blocks)
                             sbpr =[me.available_blocks.Has_subprotocol];
                             showlbl= strcat({me.available_blocks.main},{' : '},{sbpr.fulltext});
                         else
                             showlbl={};
                         end
                         ui = uicontrol('style','popupmenu','units','normalized','position',[pos+[.5,-w*.9],.4,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.block_callback(a,b),'string',[{'Blocks:'},showlbl],'HorizontalAlignment','left','value',val);
                         pos = pos+[0 -3]*w;
                    case 'Contacts'
                        if isempty(me.available_contacts) && ~isempty(me.subjid)
                            me.subject = me.subjid;
                            
                        end
                        if ~isempty(me.contacts)
                             val = find(ismember([me.available_contacts.Contact_Number],me.contacts));
                        else
                            val=1;
                        end
                        if ~isempty(me.available_contacts)
                            hlbl = [me.available_contacts.Has_Label];
                            hlbl={hlbl.fulltext};
                        else
                            hlbl={};
                        end
                          uicontrol('style','text','units','normalized','position',[pos + [.6,+1*w], .2,w*.9],'string',fields{k},'HorizontalAlignment','center');
                         ui = uicontrol('style','listbox','units','normalized','position',[pos+[.575,-5*w],.375,w*6],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.contacts_callback(a,b),'string',hlbl,'HorizontalAlignment','left','value',val,'max',2);
%                          pos = pos+[0 -4]*w;
                    case 'Response Type'
                         uicontrol('style','text','units','normalized','position',[pos+[0 -.5*w], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
                         ui = uicontrol('style','listbox','units','normalized','position',[pos+[.2,-2.5*w],.35,3.5*w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)me.response_callback(a,b),'string',{'Response Types:',me.available_response_types{:},'Create new response'},'HorizontalAlignment','left');
                         pos = pos+[0 -3]*w;
                    case 'Create New Response'
                         ui = uicontrol('style','pushbutton','string',fields{k},'units','normalized','position',[pos+[.675,-2.5*w*.8],.25,.8*w*2],'Callback',@(a,b)me.newresp(a,b));
                         pos = pos+[0 -3]*w;
               
%                           pos = pos+[0 -1]*w;
                        
                    otherwise
                          uicontrol('style','text','units','normalized','position',[pos + [0,-.025], .2,w*.9],'string',fields{k},'HorizontalAlignment','left');
                         ui = uicontrol('style','edit','units','normalized','position',[pos+[.2,0],.35,w*.9],'backgroundcolor',[1 1 1],'Callback',@(a,b)callback(fields{k},a,b),'HorizontalAlignment','left');
                         pos = pos+[0 -1]*w;
                end
                        handles.( regexprep(fields{k},'\s','_'))= ui;
            end
          
            if nargin < 2
                me.handles = handles;
            end
%              if ~isempty(me.subjid)
%                  me.subject = me.subjid;
%               end
            if ~isempty(me.guih) && ishandle(me.guih) && nargin <2 & ~isempty(me.subject)
               me.fillgui ;
%                 me.block = me.blockid;
%                 me.contacts = me.contactnums;
%                 me.response_type=me.responsetypeid;
            end
%              if ~isempty(me.blockid)
%                  me.block = me.blockid;
%              end
        end
        
        function upload(me,fname,ui,filetype,description,tofile,clobber)
            if nargin <3 
                filetype={};
            end
            if nargin < 6
                tofile = '';
            end
            if nargin < 7
                clobber=false;
            end
            if nargin < 2 || isempty(fname)
               [fname,pth] = uigetfile({['*.',filetype],filetype;'*.*','all'}); 
               q = questdlg(sprintf('Upload \n\t%s\nto the wiki?',fname));
               
               if ~strcmp(q,'Yes')
                   return
               end
               if nargin < 5 || isempty(description)
                  description =  inputdlg('Please describe the file:','max',4);
                  description = description{1}';
                  description(end+1,:)=sprintf('\n');
               end
            else
                [pth,fname,ext] = fileparts(fname);
                fname = [fname,ext];
            end
            if isempty(tofile)
                tofile=fname;
            end
            me.uploadFile(fullfile(pth,fname),tofile,description(:)',sprintf('Uploaded from the %s api interface',mfilename),clobber);
          
            if ~isempty(ui)&&nargin >2 && ishandle(ui)
                set(ui,'string',fname)
                 me.mwupdate(get(ui,'Tag'));   
            end
               
           
        end
        function mwupdate(me,fldin,varargin)
           
            fldns = regexprep(setdiff(me.fields,{'Create New Response','Submit'}),'\s','_');
            propmap = cat(1,fldns,lower(fldns));
            propmap = struct(propmap{:});
            propmap.Summary_Media='mediafile';
            propmap.Summary_PDF='pdf';
            propmap.PDF_Page='pdfpage';
            
            fld = regexprep(fldin,'\s','_');
            switch get(me.handles.(fld),'style')
                case {'popupmenu','listbox'}
                     str = get(me.handles.(fld),'string');
                     val = get(me.handles.(fld),'Value');
                     val(val<=1)=[];          
                     if isempty(val)
                        me.(propmap.(fld))='';
                     else
                        me.(propmap.(fld))=str{val};
                     end
                otherwise
                    me.(propmap.(fld))=get(me.handles.(fld),'string');
            end
            
            
           
            
        end
        function fillgui(me)
           
            fldns = regexprep(setdiff(me.fields,{'Create New Response','Submit'}),'\s','_');
            propmap = cat(1,fldns,lower(fldns));
            propmap = struct(propmap{:});
            propmap.Summary_Media='mediafile';
            propmap.Summary_PDF='pdf';
            propmap.PDF_Page='pdfpage';
            
            for k = 1:length(fldns)
                fld = regexprep(fldns{k},'\s','_');

                switch get(me.handles.(fld),'style')
                    case {'popupmenu','listbox'}
                        switch fld
                            case 'Block'
                                cmpstr = {me.available_blocks.main};
                                 val = find(ismember(cmpstr,me.block))+1;
                            case 'Contact'
                                cmpstr = [me.available_contacts.Contact_Number];
                                 val = find(ismember(cmpstr,me.contacts));
                            otherwise
                                 str = get(me.handles.(fld),'string');
                                 val = find(strcmp(me.(propmap.(fld)),str));
                        end
                         if isempty(val)
                             val=1;
                         end
                         set(me.handles.(fld),'value',val);
                    otherwise
                         set(me.handles.(fld),'string',me.(propmap.(fld)));
                end

            end 
           
            
        end
        function response_callback(me,ui,varargin)
           val = get(ui,'value')-1;
           val(val<=0)=[];
           if val > length(me.available_response_types)
              me.newresp;                
           else
              me.mwupdate('Response Type'); 
           end
               
        end
        function subject_callback(me,a,b)
            val = get(a,'Value');
            val(val==1)=[];
            me.subject = me.available_subjects{val-1};
           
        end
        function set.response_type(me,a)
           me.responsetypeid=a;
           val = find(strcmp(me.available_response_types,a))+1;
           if ~isempty(a)
               if isempty(val)
                  inp = input(sprintf('%s does not exist. Create it (Y/N)?',a),'s');
                  if strcmpi(inp,'Y')
                      me.newresp;
                  end
               end
               if ~isempty(me.guih)&&ishandle(me.guih) 
                  set(me.handles.Response_Type,'Value',val);
               end
           end
        end
        
        function a=get.response_type(me)
           a=me.responsetypeid; 
        end
        function set.subject(me,a)
             me.subjid=a;
             
            fprintf('\nLoading Subjects...');
            sresp = me.askargs({'Category','Subject','Is a fake subject','false'},{'Has ID'},{'order','descending'});
            me.available_subjects=[sresp.Has_ID];
            
            fprintf('\nLoading blocks...');            
            resp = me.askargs({'Category','Recording Block','Subject',me.subjid},{'Block','Has subprotocol'});
            sbp = [resp.Has_subprotocol];
            me.available_blocks = resp;%{resp.main};
            
            fprintf('\nLoading Response types...');

            rtresp = me.askargs({'Category','Response Type'});
            

            fprintf('\nLoading Contacts...');

            elcresp = me.askargs({'Category','Electrode Contact','Subject',me.subjid},{'Contact Number','Has Label'},{'sort','Contact Number'});
            
%             elc = [resp.Contact_Number];
            me.available_contacts = elcresp;
            me.available_response_types={rtresp.main};
            if ~isempty(me.guih) && ishandle(me.guih) && ishandle(me.handles.Subject)
                
                val = find(strcmp(me.available_subjects,me.subjid));
                
                set(me.handles.Subject,'value',val+1);
                
                
                disp = strcat({resp.main},{' : '},{sbp.fulltext});
                
                set(me.handles.Block,'string',disp);
                me.block = me.blockid;
                if isempty(get(me.handles.Block,'value'));
                    set(me.handles.Block,'Value',1)
                end
                hlbl = [elcresp.Has_Label];
                set(me.handles.Contacts,'string',{'Contacts',hlbl.fulltext});
                me.contacts = me.contactnums;
                if isempty(me.contacts)
                    cval = 1;
                else
                    cval = find(ismember([me.available_contacts.Contact_Number],me.contacts))+1;
                end
                set(me.handles.Contacts,'Value',cval);
                
            end            
        end
        function a = get.subject(me)
           a = me.subjid; 
        end
        function set.block(me,a)
             me.blockid=a;
             if ~isempty(me.guih)&&ishandle(me.guih)
                 val = find(strcmp({me.available_blocks.main},a))+1;
                 
                set(me.handles.Block,'Value',val)
             end
             
        end
        function a = get.block(me)
           a = me.blockid; 
        end
        function contacts_callback(me,a,b)
           val = get(a,'Value');
           val(val==1)=[];
           if ~isempty(me.available_contacts)
               me.contacts = [me.available_contacts(val-1).Contact_Number];
           end
        end
        function set.contacts(me,a)
             me.contactnums=a;
             if ~isempty(me.guih)&&ishandle(me.guih)
                 if ~isempty(a)
                 val = find(ismember([me.available_contacts.Contact_Number],a));
                 if isempty(val)
                     val=1;
                 end
                 
                    set(me.handles.Contacts,'Value',val+1)
                 else
                    set(me.handles.Contacts,'Value',1)
                 end    
             end
             
        end
        function a = get.contacts(me)
           a = me.contactnums; 
        end
        
        function block_callback(me,a,b)
            val = get(a,'Value')-1;
            if val==0
                me.block='';
            else
                me.block ={me.available_blocks(val).main};
            end
        end
        
        function newresp(me,varargin)
            
            fields={'Label','Stimulus','Modalities','Key Words','Description','Submit'}; %#ok<*PROPLC>
            [fig,handles] = me.gui(fields,@(varargin)[]); %#ok<*PROP>
%             dummyfig=figure('visible','off','position',[0 0 0 0]);
            set(handles.Submit,'callback',@(varargin)me.newresp_callback(fig,handles,fields,varargin));
            set(handles.Label,'string',me.response_type)
            uiwait(fig);
%             uiwait(dummyfig)
        end
        function newresp_callback(me,fig,handles,fields,varargin)
                
            fieldns = regexprep(fields,'\s*','_');
            vals = cat(1,fieldns,cellfun(@(x)get(handles.(x),'string'),fieldns,'uniformoutput',false));
            str = struct(vals{:});
            pagename = str.Label;
            txt = me.readPage(pagename);
            if ~isempty(txt)
               warning('This response type already exists!')
               rts = get(me.handles.Result_Type,'string');
               set(me.handles.Result_Type,'value',find(strcmp(str.Label,rts)))
               return
            end
            descr = str.Description';
            descr(end+1,:)=sprintf('\n');
            txt = sprintf('{{Response Type\n|Stimuli=%s\n|Modalities=%s\n|Key Words=%s\n|Description=%s\n}}',str.Stimulus,str.Modalities,str.Key_Words,deblank(descr(:)'));
            
            
            me.writeToPage(pagename,txt);
               
            delete(fig)
             me.query_params=struct('action','askargs','conditions','Category:Response Type');
              resp = me.request;
              if ~isempty(resp.query.results)
                resp = struct2cell(resp.query.results);
                resp = [resp{:}];
                me.available_response_types = {resp.fulltext};
              end
              
              if ~isempty(me.guih) && ishanlde(me.guih)
               set(me.handles.Result_Type,'value',find(strcmp(str.Label,me.available_response_types)));
              end         
             
            
        end
        function clearfields(me)
            
              fldns = regexprep(setdiff(me.fields,{'Create New Response','Submit'}),'\s','_');
            propmap = cat(1,fldns,lower(fldns));
            propmap = struct(propmap{:});
            propmap.Summary_Media='mediafile';
            propmap.Summary_PDF='pdf';
            propmap.PDF_Page='pdfpage';
            propmap.Subject='subjid';
            for k = 1:length(fldns)
                 me.(propmap.(fldns{k}))=[];
            end
        end
            
        function write(me,clobber)
            
            if ~isempty(me.page_name)
                txt = me.readPage(me.page_name);
            else
                txt = '';
            end
            if isempty(txt) || clobber
                ctxt= [sprintf('%i,',me.contacts)];
                ctxt(end)='';
                btxt= sprintf('%s,',me.block);
                btxt(end)='';
                txt = sprintf('{{Result\n|Subject=%s\n|Contacts=%s\n|Block=%s\n}}\n',me.subject,ctxt,btxt);
            else
                
                me.query_params = struct('action','browsebysubject','subject',me.page_name);
                resp = me.request;
                q = struct2cell(resp.query.data);
                q(1,:) = regexprep(q(1,:),'^_','');
                params= struct(q{:});

                if ~isequal(me.subject,params.Has_subject.item(1:3))
                    error('This page already exists with a different subject number')
                end
            end
                
%             vtxt = sprintf('%g',me.values);
            
            fldns = regexprep(setdiff(me.fields,{'Create New Response','Submit'}),'\s','_');
            propmap = cat(1,fldns,lower(fldns));
            propmap = struct(propmap{:});
            propmap.Summary_Media='mediafile';
            propmap.Summary_PDF='pdf';
            propmap.PDF_Page='pdfpage';
            fldtxt = '';
            for k = 1:length(fldns)
                val = me.(propmap.(fldns{k}));
                if isnumeric(val)
                    val=sprintf('%g,',val);
                    val(end)='';
                end
                fldtxt = sprintf('%s\n|%s=%s',fldtxt,fldns{k},val);
            end
            rtxt = sprintf('\n{{Response%s\n}}',fldtxt);
            wrtxt = sprintf('%s%s',txt,rtxt);
            
            if isempty(me.page_name)                
                me.page_name = sprintf('%s/Contact %s Results',btxt,ctxt);
            end
%             me.writeToPage(me.page_name,wrtxt,1-clobber)
            me.writeToPage(me.page_name,wrtxt,0)
            
        end
        
        function txt = read(me)
           
            
            txt = me.readPage(me.page_name);
            if ~isempty(txt)  && ~me.overwrite
                params= me.query_params;
                
                me.subject = params.Has_subject;
                me.block = unique({me.block,params.Has_block});
                me.contacts = unique({me.contacts,params.Has_electrode_contact});
                
            end
            
        end
        function q = query_page(me)
            
            
            me.query_params = struct('action','browsebysubject','subject',me.page_name);
            res = me.request;
            if isfield(res,'query')
                q = struct2cell(res.query.data);
                dd = arrayfun(@(x){regexprep({x.dataitem.item},'#\d*#','')},res.query.data,'uniformoutput',false);
                q(1,:) = regexprep(q(1,:),'^_','');
                q(2,:) = dd;
                q = struct(q{:});
            else
                q=[];
            end
        end    
            
            
    end
        
        
        
end
   