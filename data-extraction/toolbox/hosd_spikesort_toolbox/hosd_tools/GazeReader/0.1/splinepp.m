function [R,Cmat,Yout] = splinepp(Y,T,nknots,tt,varargin)

% [R,Cmat,Yout] = splinepp(Y,T,nknots,tt)
%
% Create a spline model based on the distribution of events within some
% window of obervation. This employs the format used by mnlfit: the null 
% outcome is modeled explicitly in the design matrix. Yout indicates the
% row associated with the outcome.
%
% See also MNLFIT


if nargin < 3 || isempty(nknots)
    nknots = 4;
end

if nargin < 4 || isempty(tt)
    tt =(0:size(T,1)-1)';
end
%%
dY = [0;diff(Y)]; % Detect edges in case Y is based on simple thresholding

N = cumsum(sum(abs(dY(T)) + 1./size(T,1),2)); 

N = round(N./repmat(N(end,:),size(N,1),1)*nknots);

knots = [tt(1)-diff(tt(1:2))*.5;tt(unique(find(diff(N))));tt(end)];

T = T(:,randperm(size(T,2))); %The order of columns in T is randomized here in case the windows overlap. This randomizes whether a given point is
                              %assigned to the earlier or later window.
x = zeros(size(Y));
x(T(:)) = repmat(tt,size(T,2),1);
xintcpt = zeros(size(x));
xintcpt(T(:))=1;

[R,Cmat] = buildsplinereg(kron(x,[0 1]'),knots,'noptions',2,'codeincr',1,'postmultiply',kron(xintcpt,[0 1]'),varargin{:});
% R(2) = makeregressor(kron(xintcpt,[0 1]'),'codeincr',2,'label','intercept','noptions',2);
% Cmat(:,end+1)=0;

Yout = zeros(2,length(Y));
Yout(1+Y+2*(0:length(Y)-1)')=1;


