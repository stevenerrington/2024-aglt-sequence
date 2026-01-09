function [chanceThr, nullAcc] = svm_decode_permutation(X, Y, nFolds, nPerm)
% Permutation-based chance threshold for SVM decoding
%
% Outputs:
%   chanceThr : 95th percentile of null distribution
%   nullAcc   : null accuracy distribution

nullAcc = nan(nPerm,1);

for p = 1:nPerm
    Yshuff = Y(randperm(numel(Y)));
    nullAcc(p) = svm_decode_cv(X, Yshuff, nFolds);
end

chanceThr = prctile(nullAcc, 95);

end
