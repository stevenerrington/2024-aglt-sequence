function [acc, confMat, classAcc] = svm_decode_cv(X, Y, nFolds)
% Linear SVM with cross-validation
% Returns:
%   acc      - overall accuracy
%   confMat  - confusion matrix (rows = true, cols = predicted)
%   classAcc - per-class accuracy

cv = cvpartition(Y,'KFold',nFolds);

allTrue = [];
allPred = [];

for f = 1:cv.NumTestSets
    trainIdx = training(cv,f);
    testIdx  = test(cv,f);

    model = fitcecoc( ...
        X(trainIdx,:), Y(trainIdx), ...
        'Learners','linear', ...
        'Coding','onevsall');

    yhat = predict(model, X(testIdx,:));

    allTrue = [allTrue; Y(testIdx)];
    allPred = [allPred; yhat];
end

% Overall accuracy
acc = mean(allPred == allTrue);

% Confusion matrix
confMat = confusionmat(allTrue, allPred);

% Per-class accuracy (diagonal / row sum)
classAcc = diag(confMat) ./ sum(confMat,2);

end
