
function [out,bl] = rmbaseline(bx,usepts)


% Fits a polynomial to the spectrum and normalizes by the fitted value.

% C Kovach 2013
% 
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/DBT/rmbaseline.m $
% $Revision: 491 $
% $Date: 2014-05-06 17:57:17 -0500 (Tue, 06 May 2014) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 2 || isempty(usepts)
    usepts = true;
end
polyord = 8;

mn = mean(abs(bx.blrep),1);
kp = find(mn~=0 & usepts);
p = polyfit(kp./length(mn),log(mn(kp)),polyord);

bl = polyval(p,(1:length(mn))/length(mn));



out = bx.blrep./repmat(exp(bl),size(bx.blrep,1),1);

