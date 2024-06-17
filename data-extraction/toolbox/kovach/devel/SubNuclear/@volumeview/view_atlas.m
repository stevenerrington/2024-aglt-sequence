% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/@volumeview/view_atlas.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

function view_atlas
persistent atlas AP tmpl2maimat 


if isempty(atlas)
    load maiatlas2 atlas
    load mai_template_mni_aligned

end


