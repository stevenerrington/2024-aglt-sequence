function svn_add(fnames)

%
%
%    svn_add_all_files( files )
%
%  Adds specified files to the repository.
%
%  For m-files revision data is also inserted after the first block of
%  comments.
%
% See also insert_svn_revision_data Add_this_to_matlab_scripts_for_revision_data.readme
%
%


%  C Kovach 2011
%
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/svntools/svn_add.m $     
% $Revision: 245 $
% $Date: 2013-07-23 17:17:37 -0500 (Tue, 23 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------
%

if ~iscell(fnames)
     fnames = {fnames};
end
for i = 1:length(fnames)
    
    [a,b,c] = fileparts(fnames{i});

    

    com = sprintf('svn add %s',fnames{i});

    system(com);

    if isequal(c,'.m')
        svn_insert_revision_data;
    end

end



