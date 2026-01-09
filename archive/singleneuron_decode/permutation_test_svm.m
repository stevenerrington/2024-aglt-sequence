function pVal = permutation_test_svm(X, Y, nFolds, accReal, nPerm)
% Compute permutation-based p-value
nTrials = numel(Y);
permAcc = nan(nPerm,1);

for k = 1:nPerm
    Yperm = Y(randperm(nTrials));   % shuffle labels
    permAcc(k) = svm_decode_cv(X, Yperm, nFolds); % reuse existing function
end

% p-value = fraction of permutations >= real accuracy
pVal = mean(permAcc >= accReal);
end