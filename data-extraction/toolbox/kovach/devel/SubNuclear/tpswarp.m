
function [warpfun,d] = tpswarp(X1,X2,reg)

% Warpfun = tpswarp(Xfrom,Xto, reg)
% Thin-plate spline warping. 
%      
%      Xfrom : row matrix of coordinates in the starting space.
%      Xto:   coordinates in the target space.
%      reg:  regularization term (default reg=0). For large values warping
%            approaches an affine transformation and for reg=0 its fits
%            the control points exactly.
%
%      Warpfun is a function handle where
%      Xto == Warpfun(Xfrom)   ( for reg = 0 );
% 
% [Warpfun,A] = tpswarp(Xfrom,Xto, reg)
%      Also returns the affine component of the transformation in matrix A.

% C Kovach 2013

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/tpswarp.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


if nargin < 3 || isempty(reg)
    reg = 0;
end

squash = @(x)sqrt(sum(cat(3,x{:}),3));
dxc = @(x,c,k) (repmat(x(:,k),1,length(c(:,k)))-repmat(c(:,k)',size(x,1),1)).^2;
mdxc = @(x,c) squash(arrayfun(@(k)dxc(x,c,k),1:size(x,2)-1,'uniformoutput',false));
tpsmat = @(x,c) mdxc(x,c).^2.*log(mdxc(x,c)+eps);
pad1 = @(x) cat(2,x,ones(size(x(:,1))));
KX = @(x)tpsmat(pad1(x),X1);

[Q,R] = qr(pad1(X1));
Q1 = Q(:,1:size(X1,2)+1);
Q2 = Q(:,size(X1,2)+2:end);
R = R(1:size(X1,2)+1,:);
TPScoef = Q2*(((Q2'*KX(X1)*Q2+reg*eye(size(Q2,2))))'\(Q2'*X2));
d = R^-1*Q1'*(X2-KX(X1)*TPScoef);
warpfun = @(x) pad1(x)*d + KX(x)*TPScoef;

