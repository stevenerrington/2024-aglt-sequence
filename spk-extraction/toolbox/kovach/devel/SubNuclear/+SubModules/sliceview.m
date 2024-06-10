
classdef sliceview < SubModules.module

%
% The slice module generates a planar view of arbitrary orientation. 
%
% To use sliceview, select 3 or more points or 1 or more meshes containing
% at least 3 vertices in the respective menus, then double click on the 
% 'sliceview' module. A new figure window and axis will appear.
%
% The orientation of view is determined from the first two principle components 
% of the group of selected points. 
%
% An extra line will appear in the crosshars for all other open axes,
% indicating the lines of intersection for the plane.
%
% You can create a view of arbitrary orientation by selecting 3 points in 
% the desired plain. 
%
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/+SubModules/sliceview.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2013

    methods

        function me = sliceview(vv,X,addsis,varargin)        
            me.initialize(vv,varargin{:});
            vv = me.parent;
            cro = vv.current_object;
            fld = class(cro);
            indx = ismember([vv.(fld)],cro);

            if nargin < 2 || isempty(X)
                
                switch class(cro)
                    case {'points','lines'}


                        X = cat(1,vv.points(indx).coord);

                    case 'meshes'
                        decell = @(x)cat(1,x{:});
                        X = decell( arrayfun(@(x) x.trirep.X, vv.meshes(indx),'uniformoutput',0));
%                     case {'volumes','transforms'}
                        
                    otherwise
                        warning('To use %s first select points or meshes used to define the plane of view.',mfilename)
%                         warning('Showslice does nothing for %s',class(cro))
                        return          

                end
            end
                
            if nargin < 3 
                addsis = true;
            end

            if size(X,1) < 3
                error('Requires at least 3 points');
            end

            pt = mean(X);
            vv.current_point = pt;
            [u,~]= svd(cov(X));
            nv = u(:,3)';
            vc = size(vv.current.volumes.image.Data)/2;
             p0 = ((pt-vc)*nv')*nv+vc;
            [sl,axv,trmat] = vv.makeslice(nv,p0);
            nfig = figure;
            vv.addfig(nfig);
            axh = axes;
            ax=vv.addaxis(axh); 
            szsl = size(sl);
            szsl(end+1) = 0;
            ax.axdim = szsl([3 2 3 1]);

            % axi = find(ax==vv.axes);

            ax.normvec=nv;
%             ax.axvec=axv;
            ax.basetrmat=trmat;
            vv.resetaxis(ax); 
%             if addsis
%                 for i = 1:length(vv.sisters)
%                     SubModules.sliceview(vv.sisters(i),X,false);
%                 end
%             end
            me.data.ax = ax;
            vv.plotupdate();
        end
        
        function update(me)
           % Updating already handled 
        end
            
        
    end
end