% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/SubNuclear/@volumeview/view_atlas.m $
% $Revision: 258 $
% $Date: 2013-07-25 10:38:04 -0500 (Thu, 25 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------

function view_atlas
persistent atlas AP tmpl2maimat 


if isempty(atlas)
    load maiatlas2 atlas
    load mai_template_mni_aligned

end


