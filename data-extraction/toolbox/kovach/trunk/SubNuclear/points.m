% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/SubNuclear/points.m $
% $Revision: 1153 $
% $Date: 2019-02-03 21:25:17 -0600 (Sun, 03 Feb 2019) $
% $Author: ckovach $
% ------------------------------------------------

classdef points < plotobj
    properties
            coord = [];
            file;
            path;
            units = 'vox';
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
            if nargin > 0 && ~isempty(varargin)>0 && isa(varargin{1},'points')
                fln = fieldnames(varargin{1});
                
                for i = 1:length(fln)
                    for ii = 1:length(varargin{1});
                        imo(ii).(fln{i}) = varargin{1}(ii).(fln{i}); 
                    end
                end
            else
                for i = 1:2:length(varargin)
                    imo.(varargin{i})=varargin{i+1};
                end
            end
        end
    end
end
