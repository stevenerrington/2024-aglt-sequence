function overlap = compute_subspace_overlap_cv(X1, X2, nPCs, nFolds)
% Cross-validated subspace overlap between X1 and X2
% Inputs:
%   X1, X2 : [neurons x trials/timepoints] matrices
%   nPCs   : number of PCs to use
%   nFolds : number of cross-validation folds
% Output:
%   overlap : scalar between 0 and 1

% Default number of folds
if nargin < 4
    nFolds = 5;
end

% Ensure neurons x samples orientation
if size(X1,1) < size(X1,2), X1 = X1'; end
if size(X2,1) < size(X2,2), X2 = X2'; end

% Subtract mean
X1 = bsxfun(@minus, X1, mean(X1,2));
X2 = bsxfun(@minus, X2, mean(X2,2));

% Partition indices
nSamples = size(X1,2);
indices = crossvalind('Kfold', nSamples, nFolds);

varExplained = zeros(nFolds,1);

for i = 1:nFolds
    % Split train/test
    testIdx = (indices == i);
    trainIdx = ~testIdx;

    % PCA on X1 training set
    [U1, ~, ~] = svd(X1(:, trainIdx), 'econ');
    basis1 = U1(:, 1:nPCs);

    % Project test X2 data onto basis
    proj = basis1' * X2(:, testIdx);

    % Compute variance explained
    var_proj = sum(var(proj, 0, 2));
    var_X2   = sum(var(X2(:, testIdx), 0, 2));

    varExplained(i) = var_proj / var_X2;
end

% Average across folds
overlap = mean(varExplained);
end
