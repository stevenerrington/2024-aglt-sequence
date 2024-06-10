% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/contactplot.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

function pl = contactplot(COR,C,varargin)

hold on
for i = 1:size(COR,1)
    
    pl(i) = plot3(COR(i,1),COR(i,2),COR(i,3),varargin{:});
    
    set(pl(i),'color',cm2rgb(C(i)));
end


