% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/lines.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

classdef lines < points
     methods
        function imo = points(varargin)
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end
        end
    end
end
