
classdef module < handle  
    
% Abstract class definition for a SubNuclear module.
%
% To create your own module add a file [your_module].m to the +module
% directory according to the following template:
%
%  classdef your_module < SubModules.module
%
%     methods
%        function this = your_module(vv)  % constructor
%           
%           this.parent = vv;  % The callling volumeview object must be contained in
%                              % this.vv
%          
%           ... (your code)
%
%         end
%
%       function update(me)
%
%           ... (Code to run, if any, after every plot update. If there is
%                nothing to do leave this empty, but still included it)     
%
%       end
%
%    end
% end
%
% See also Classdef
    
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/transforms.m $
% $Revision: 279 $
% $Date: 2013-07-26 12:29:39 -0500 (Fri, 26 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------

    properties
        label = '';
        notes= '';
        objectid = -1;
        parent;
        data
    end
    
    methods
        function initialize(me,varargin)
            if nargin < 2 || isempty(varargin{1})
                vv = volumeview([]);
            else
                vv = varargin{1};
            end
            me.parent = vv;
            if ~isempty(vv)
                me.objectid = vv.objectiter;
            end
        end
        function updateview(me)
            for i = 1:length(me)
                me(i).update();
            end
        end

    end
    methods (Abstract)     
        update(me);
    end
  
end