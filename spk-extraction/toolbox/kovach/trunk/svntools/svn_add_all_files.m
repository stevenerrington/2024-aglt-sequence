function svn_add_all_files(ext)

%
%
%    svn_add_all_files( ext )
%
%  Adds all files of type ext (default '.m') to the repository.
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
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/svntools/svn_add_all_files.m $     
% $Revision: 255 $
% $Date: 2013-07-23 22:12:22 -0500 (Tue, 23 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------
%

if nargin < 1 || isempty(ext)
     ext = 'm';
end

[a,b,c] = fileparts(ext);

if ~isempty(a)
    ext = fullfile(a,[b,c]);
elseif isempty(c)
    ext = ['*.',b];
    c = ['.',b];
end


com = sprintf('svn add %s',ext);

system(com);

if isequal(c,'.m')
    svn_insert_revision_data;
end




