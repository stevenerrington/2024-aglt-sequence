% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/SubNuclear/meshes.m $
% $Revision: 891 $
% $Date: 2017-06-29 09:57:38 -0500 (Thu, 29 Jun 2017) $
% $Author: ckovach $
% ------------------------------------------------

classdef meshes < plotobj
    properties
            trirep = [];
            file;
            path;
    end
    
     methods
        function imo = meshes(varargin)
            if nargin> 0 && isa(varargin{1},'meshes')
                fldn = fieldnames(varargin{1});
                fldn = setdiff(fldn,{'ploth','objectid'});
                for k = 1:length(fldn)
                    imo.(fldn{k}) = varargin{1}.(fldn{k});
                end
            elseif nargin>0
                
                for i = 1:2:length(varargin)
                    imo.(varargin{i})=varargin{i+1};
                end
            end
        end
    end
end
