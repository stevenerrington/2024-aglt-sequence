% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/projectcontact.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------


function [facets,mnd] = projectcontact(C,bem,layers)

player = ismember(bem.layer,layers);

for i = 1:size(C,1)
    
    D = sqrt(sum((bem.FC - repmat(C(i,:),size(bem.FC,1),1)).^2,2));
    D(~player) = Inf; 
    [mnd(i),facets(i)] = min(D);
end
