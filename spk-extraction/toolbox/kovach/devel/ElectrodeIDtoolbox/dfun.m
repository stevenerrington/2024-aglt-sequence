
function D = dfun(A,B)

% D = dfun(A,B)
%Compute distances between pairs of points in A and B

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/dfun.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

NB = size(B,1);
NA = size(A,1);
D = zeros(size(A,1),size(B,1));
for i = 1:size(B,2)
    D = D + (repmat(A(:,i),1,NB) -repmat(B(:,i)',NA,1)).^2;
end
D = sqrt(D);