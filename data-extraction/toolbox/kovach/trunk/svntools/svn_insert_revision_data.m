function insert_svn_revision_data(fnames)

%
%
%    insert_svn_revision_data
%
% Scans all m-files in the current path and inserts a block of commented text
% after the first block of comments, which contains SVN keywords needed to 
% automatically insert revision information. If the file contains no comments, 
% then the text is inserted at the top. The keywords are replaced with revision
% information every time the file is updated or checked out.
%
%    insert_svn_revision_data  file.m 
%
% Operate on specified file,
%
%
%    insert_svn_revision_data  path
%
% Operate on files in given path.
%
%
% insert_svn_revision_data({file1.m file2.m ...}), 
% 
%    Operate on given files
%
%
% See also Add_this_to_matlab_scripts_for_revision_data.readme
%
%


%  C Kovach 2011
%
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/svntools/svn_insert_revision_data.m $     
% $Revision: 350 $
% $Date: 2013-09-24 17:02:17 -0500 (Tue, 24 Sep 2013) $
% $Author: ckovach $
% ------------------------------------------------
%

text ={'% ----------- SVN REVISION INFO ------------------'   %Don't change this line
       '% #URL#'            %Dollar signs replaced with # so that svn doesn't meddle here.
       '% #Revision#'
       '% #Date#'
       '% #Author#'
       '% ------------------------------------------------'};
   
textrep = regexprep(sprintf('%s\n',text{:}),'#','$');



if nargin < 1
    
    d = dir('*.m');
    fnames = {d.name};
    
else
    
    if ischar(fnames)  % If input is not a cell array of file names...
        [p,f,e] = fileparts(fnames); %#ok<*NASGU>
        if ~isempty(e)  % assume input is file name
%             if isempty(e)
%                 e='.m';
%             end
            fnames = {fullfile(p,[f,e])}; 
        else   %  assume input is a path.
            p=fullfile(p,f);
            if strncmp(p,'*',1);
                d = dir(p);
                fnames = {};
                for i = 3:length(d)
                    if d(i).isdir && ~strncmp(d(i).name,'.',1)
                        insert_svn_revision_data(fullfile(p,d(i).name,filesep));
                    end
                end
            else
                d = dir(fullfile(p,'*.m'));
                fnames = strcat(p,filesep,{d.name});
            end
        end
    end    
    
end
    
for i = 1:length(fnames)
    
    com = sprintf('svn propset svn:keywords "Author Date Revision URL Id" %s',fnames{i});
    [status,result] = system(com); %#ok<*ASGLU>
    if ~isempty(regexp(result,'svn: warning:','once'));   
        continue
    end
    
    fid = fopen(fnames{i},'r');
    fcont = fread(fid,Inf,'uchar=>char')';
    fclose(fid);
    vin = regexp(fcont,'- SVN REVISION INFO -','once');
    if ~isempty(vin) % file already contains info
        continue
    end
    
    commentblk = regexp(fcont,'[ \t]*%[^\n]*\n[ \t]*[^%\n]*\n','once','end'); % find first end of a comment
    
    if isempty(commentblk)
        newfcont = sprintf('%s\n%s',textrep,fcont);
    else
        newfcont = sprintf('%s%s\n%s',fcont(1:commentblk),textrep,fcont(commentblk+1:end));
    end        
    
    %Now write the updated file
    fid = fopen(fnames{i},'w');
    fwrite(fid, newfcont);
    fclose(fid);
    
end



