function inds = randpick(A,n)
    
%Randomly pick indices in A based on the relative magnitudes of the values
%of A.
if nargin < 2
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

    n = 1;
end
if min(size(A)) == 1
 A = A(:);
 A = repmat(A,1,n);
end

csA = cumsum(A)./repmat(sum(A),size(A,1),1);

randvals = rand(1,size(A,2));

diffs = csA - randvals(ones(size(csA,1),1),:);
diffs(diffs < 0) = Inf;

[mn,inds] = min(diffs);