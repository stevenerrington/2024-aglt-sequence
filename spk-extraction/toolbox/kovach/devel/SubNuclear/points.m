% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/points.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

classdef points < plotobj
    properties
            coord = [];
            file;
            path;
    end
     methods
        function imo = points(varargin)
            if nargin > 0 && isnumeric(varargin{1})            
                X = varargin{1};
                varargin(1) = [];
                if size(X,1) == 1
                    imo.coord = X;
                else
                    for i = 1:size(X,1)
                        imo(i) = points(X(i,:),varargin{:});
                    end
                end
            end
            for i = 1:2:length(varargin)
                imo.(varargin{i})=varargin{i+1};
            end

        end
    end
end
