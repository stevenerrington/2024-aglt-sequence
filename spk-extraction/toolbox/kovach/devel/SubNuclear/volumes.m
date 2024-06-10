% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/volumes.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

classdef volumes < plotobj
    properties
           image = medimage([]);
           intensity_range = [0 255];
           file;
           path;
    end
    
     methods
        function imo = volumes(varargin)
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end
        end
        function a = get.file(me)
            a = me.image.file;
        end
        function a = get.path(me)
            a = me.image.path;
        end
    end
end
