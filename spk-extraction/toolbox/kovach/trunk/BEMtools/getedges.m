
function bem = getedges(bem)

%Adds a Nedge x 2 matrix of edges and Nface x 1 vector which maps into the 
%edges matrix toe  the bem structure. 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/getedges.m $
% $Revision: 37 $
% $Date: 2011-06-04 23:17:29 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------


alledges = sort(cat(1,bem.tri(:,1:2),bem.tri(:,2:3),bem.tri(:,[1 3])),2);

edges = unique(alledges,'rows');

tridim = [1 2;2 3; 3 1];

for i = 1:3    
    [ism,edgemap(:,i)] = ismember(sort(bem.tri(:,tridim(i,:)),2),edges,'rows');
end

bem.edges = edges;
bem.edgemap = edgemap;