function C = ftconv(X1,X2,correl,center)

%  C = ftconv(X1,X2,correl,center)
%
%Simple, fast, and infinitely useful convolution or correlation of the columns 
% of two matrices by multiplication of their respective fourier transforms.
%
% X1 - matrix 1
% X2 - matrix 2
%
% correl = true for correlation. Default is false.
%
% center = [c1 c2]: c1 - whether to circular shift the first matrix so that the 
% middle is centered on the first element, c2 - whether to center the
% second matrix likewise.   Default is [false true].
%
% Note that this algorithm applies convolution in the native fourier space,
% which is circular. This means that the beginning of the signal will be smeared
% into the end. If this is a problem for you, then remember to zero pad the
% columns of the matrices by the size of the smaller of the two matrices.
%

% C Kovach 2010



% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/signaltools/ftconv.m $
% $Revision: 47 $
% $Date: 2011-06-30 13:07:10 -0500 (Thu, 30 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 3
    correl = false;
end
if nargin < 4
    center = [0 1];
end

sz1 = size(X1);
sz2 = size(X2);
csize = max(cat(1,sz1,sz2));

for i = 1:length(sz1)
    ind1{i} = 1:sz1(i);
    ind2{i} = 1:sz2(i);
end

A1 = zeros(csize);
A1(ind1{:}) = X1;
if center(1) 
    A1 = circshift(A1,-floor(sz1./2));
end

A2 = zeros(csize);
A2(ind2{:}) = X2;
if center(2)
        A2 = circshift(A2,-floor(sz2./2));
end


if ~correl
    C = ifftn(fftn(A1).*fftn(A2));
else
    C = ifftn(fftn(A1).*conj(fftn(A2)));
end    




