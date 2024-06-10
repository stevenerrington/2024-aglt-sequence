
function labarchive2labwiki(lapi,notebook,trees,update,overwrite)

timestamp_protocol_bug = true; %Timestamp and protocol attributes are swapped;
if nargin < 5 || isempty(overwrite)
    overwrite = false; %Overwrite entire entry 
end
if nargin < 4 || isempty(update)
    update = true; %Update any changes more recent than the last revision date (note, newer changes will not be updated if the page was modified in any way)
end
if nargin < 2
    notebook = lapi.notebooks(1);
end
if nargin < 3
    gettree = find(arrayfun(@(x)~isempty(regexp(x.label,'\d\d\d')),notebook.tree));
    trees = notebook.tree(gettree);
end

fix_ampersand = true;

tag_new_with='lwapi-labarchive_dump_new';
tag_mod_with='lwapi-labarchive_dump_edit';

persistent mwapi
if isempty(mwapi)
    mwapi = mwapi2;
end
if ~mwapi.logged_in
   mwapi.login; 
end

for sbi = 1:length(trees)

    tree = trees(sbi);
    dd = lapi.get_tree(notebook,Inf,tree);
    dats = dd(strmatch('Block',{dd.label}));
    
    if isempty(dats)
        dats = struct('branch',dd);
    end
  
    try
        sid =  regexp(dats.branch(end).entry(1).parsed.subjectnumbermandatory,'\d\d\d','match','once');
        
        mwblks = mwapi.askargs({'Category','Recording Block','Subject',  sid});
        mwblks = {mwblks.main};
    catch
        mwblks={};
    end
    
    for k = 1:length(dats.branch)
        %%
        entries = dats.branch(k).entry;
        if ~isfield(entries,'parsed') || isempty([entries.parsed])
            continue
        end
        d = entries(1).parsed;
        sid = d.subjectnumbermandatory;
        sid = regexp(sid,'\d\d\d','match','once');
        if isfield(d,'blocvkid') && ~(update || overwrite) && any(strcmp(d.blocvkid,mwblks)) 
            continue
        end
        %%
        if isempty(d.recordingsystem)
            d.recordingsystem = 'Neuralynx';
        end
        if ~isfield(d,'date')
           d.date = d.exp_date; 
        end
        if isempty(d.date)
            d.date = d.datedate;
        end
        if ~isempty(d.date)
            try
                dnum = datenum(d.date);
                d.date = datestr(dnum,'mm/dd/yyyy');
            catch
                d.date = d.datedate;
            end
        end

        if ~isfield(d,'blocvkid')
            d.blocvkid=d.blockID;
            if timestamp_protocol_bug
                proto=d.timestamp;
                d.timestamp=d.protocol;
                d.protocol=proto;
            end
        end
        blid = regexp(d.blocvkid,'-(.*)','tokens','once');
        if isempty(blid)
            blid = {sprintf('%03i',str2double(d.blocvkid))};
            d.blocvkid = sprintf('%s-%03i',sid,str2double(d.blocvkid));
        
            
        end
        block = d.blocvkid;
        if ~isfield(d,'memobox')
            if ~isfield(d,'ExperimentDetails')
                d.ExperimentDetails=d.experiment_details;
            end
           d.memobox = sprintf('%s\n%s\n%s',d.ExperimentDetails,d.RecordingNotes,d.responseNotes); 
        end
        d.memobox = regexprep(d.memobox,'\\\\[nr]','\\n');
        d.memobox = regexprep(d.memobox,'\\n',sprintf('\n'));
        
        lwstr = struct('subject',sid,'block',blid{1},'date',d.date,'time',d.timestamp,...
                   'Subprotocol',d.protocol,'Montage',sprintf('Subject/%s/Montage',sid),'System',d.recordingsystem,...
                   'Notes',d.memobox);

        lwfld = fieldnames(lwstr);
        lws = cellfun(@(x)sprintf('\n|%s=%s',x,lwstr.(x)),lwfld,'uniformoutput',false);

        page_name = d.blocvkid;
        [mw,rvdat] =mwapi.readPage(page_name);
                
        modt = datenum(cat(1,entries.modified),'yyyy-mm-ddTHH:MM:SSZ');
        [~,mxi] =max(modt);

        if ~isempty(rvdat)
             lwmodt = datenum(rvdat.timestamp,'yyyy-mm-ddTHH:MM:SSZ');
        else
            lwmodt = -Inf;
        end
        fxamp=fix_ampersand && ~isempty(regexp(d.memobox,'&','once'));
        if fxamp
            mwnotes = regexp(mw,'|\s*Notes\s*=.*?}}','match','once');
            fxamp = isempty(regexp(mw,'&','once')); %If ampersand is in labwiki notes, assume the correction was already applied
            if fxamp
                fprintf('\n%s: memobox contains an ampersand. Will overwrite entry',d.blocvkid);
            end
        end
        if isempty(mw) || overwrite || fxamp
           labwikitxt = sprintf('{{Recording Block%s\n}}',[lws{:}]);
           overwrite_this=true;
            fprintf('\nWriting %s',d.blocvkid);
               tag=tag_new_with;

        elseif  (modt(1)>lwmodt && update)
            nulabwikitxt = mw;
            for kk = 1:length(lwfld)
               nulabwikitxt = regexprep(nulabwikitxt,sprintf('\\|\\s*%s\\s*=[^\\n\\|}]*',lwfld{kk}),sprintf('|%s=%s',lwfld{kk},lwstr.(lwfld{kk})));
            end
            labwikitxt=nulabwikitxt;
            fprintf('\nUpdating %s',d.blocvkid);
   
            overwrite_this = true;
            tag=tag_mod_with;
        else
            labwikitxt = mw;
            overwrite_this = true;
            tag=tag_mod_with;
            
        end
        if ~isempty(d.timestamp) && (overwrite||isempty(regexp(mw,'\|\s*time\s*=','once')))&&isempty(regexp(labwikitxt,'\s*time\s*=\s*\d+'))
            dtstr = regexp(labwikitxt,'|\s*date\s*=([^\n]*)','tokens','once');
            dt= datenum(dtstr{1});
            if dt==floor(dt) 
               if ~isempty(regexp(dtstr{1},'00:00','once'))
                   labwikitxt = regexprep(labwikitxt,'\|\s*date\s*=\s*([^\s]*)\s+0[0:]+','|date=$1');
               end
               labwikitxt=regexprep(labwikitxt,sprintf('{{Recording Block\\s*\n'),sprintf('{{Recording Block\n|time=%s\n',d.timestamp));
            end
        end
        if ~isempty(entries(1).entryURL) && (overwrite ||isempty(regexp(mw,'\|\s*LabArchiveURL\s*=','once')))||fxamp
            %labwikitxt=regexprep(labwikitxt,sprintf('{{Recording Block\\s*\n'),sprintf('{{Recording Block\n|LabArchiveURL=%s\n',entries(1).entryURL));
            labwikitxt=regexprep(labwikitxt,'({{Recording Block.*?)\}\}',sprintf('$1|LabArchiveURL=%s\n}}',entries(1).entryURL));
        end
         if ~isempty(d.memobox) && isempty(regexp(labwikitxt,'\|\s*Notes\s*=','once'))
              labwikitxt=regexprep(labwikitxt,'({{Recording Block.*?)\}\}',sprintf('$1|Notes=%s\n}}',d.memobox));
         end
        
        re = regexp(labwikitxt,'(.*)({{\s*#set:.*?}}.*)','tokens','once');
       
        
        if ~isempty(re) && all(modt<lwmodt) && ~cellfun(@(x)isempty(x),regexp(re{2},'#set:(.*)\n}}','tokens','once'))
            labwikitxt = re{1};
             latxt = re{2};
        elseif  update || fxamp
            if isempty(re)
                latxt = sprintf('\n\n{{#set:\n}}');
            else
                 labwikitxt=re{1};
                 latxt=re{2};
            end
           for kk = 1:length(entries)    
               d = entries(kk).parsed;
               if isempty(d) || ischar(d) || modt(kk)<lwmodt && (~isempty(re) &&~cellfun(@(x)isempty(x),regexp(re{2},'#set:(.*)\n}}','tokens','once')))
                   continue
               end
                dfld = setdiff(fieldnames(d),'memobox','stable');
                
                for fldk = 1:length(dfld)
                    if strcmp(d.(dfld{fldk}),'-')
                        d.(dfld{fldk})='';
                    end
                    nuparam = sprintf('Has labarchive %s',dfld{fldk});
                    nutxt = sprintf('|Has labarchive %s=%s',dfld{fldk},d.(dfld{fldk}));
                    re2 = regexp(latxt,sprintf('\\|\\s*%s\\s*=[^\\n\\|\\}]*',nuparam),'once');
                    if ~isempty(re2)
                        if isempty(d.(dfld{fldk}))
                            latxt = regexprep(latxt,sprintf('\\n\\|\\s*%s\\s*=[^\\n\\|\\}]*',nuparam),'');
                            fprintf('\nUpdating %s: Removing %s',block,dfld{fldk});
           
                        else
                            latxt = regexprep(latxt,sprintf('\\|\\s*%s\\s*=[^\\n\\|\\}]*',nuparam),nutxt);
                            fprintf('\nUpdating %s: Replacing %s',block,dfld{fldk});
                        end
                    elseif ~isempty(d.(dfld{fldk}))
                         latxt = regexprep(latxt,'{{#set:',sprintf('{{#set:\n%s',nutxt));
                         fprintf('\nUpdating %s: Adding %s',block,dfld{fldk});

                    end
                end
                if ~isempty(mw)
                    mwprot = regexp(mw,'\|\s*Subprotocol\s*=([^\n]*)','tokens','once');

                    if isfield(d,'protocol') && ~strcmpi(deblank(mwprot{1}),deblank(d.protocol))
                        nuparam='Has labarchive protocol name conflict';
                         nutxt = sprintf('|Has labarchive protocol name conflict=%s',d.protocol);
                        re2 = regexp(latxt,sprintf('\\|\\s*%s\\s*=[^\\n\\|\\}]*',nuparam),'once');
                        if isempty(re2)                                        
                          latxt = regexprep(latxt,'{{#set:',sprintf('{{#set:\n%s',nutxt));
                          fprintf('\nWiki and Archive protocol names are different:\n\tWiki:%s\n\tArchive:%s',deblank(mwprot{1}),deblank(d.protocol));
                        end
                    end
                end
           end
        %    latxt = sprintf('{{#set: %s\n}}',[addt{:}]);
            latxt = regexprep(latxt,'\\[nr]',sprintf('\n'));
           
        else
             latxt='';
        end
        addtxt = deblank(sprintf('%s%s',labwikitxt,latxt));
%         if isempty(addtxt)
%             fprintf('%s already exists, skipping',page_name)
%             continue
%         end
         addtxt = regexprep(addtxt,'\\\\[nr]','\\n');
         addtxt = regexprep(addtxt,'\\n',sprintf('\n'));
        if ~isequal(mw,addtxt)
            comment=sprintf('Imported on %s from Labarchives entry %s (modified %s by %s)',datestr(now),char(mwapi.base64dec(entries(mxi).entryID)),entries(mxi).modified,entries(mxi).last_modified_by );
            res = mwapi.writeToPage(page_name,addtxt,1-overwrite_this,comment,tag);
        
            if isfield(res,'edit') && isfield(res.edit,'title')
                fprintf('\n...%s',res.edit.title,res.edit.result)
            else
                fprintf('\n... Error writing to %s',page_name)
            end      
            ress{k}=res; %#ok<AGROW,NASGU>
        else
            fprintf('\nNo changes in %s',block)
        end
    end
 end