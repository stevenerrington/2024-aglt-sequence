% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/plotobj.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

classdef plotobj < imobj  
    properties
        ploth = [];
        trmat = [];
        plotcolor =[];
        linestyle = '';
        marker = '';
        plotargs = {};
        show = false;
    end
end