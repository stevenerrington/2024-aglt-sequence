function acc_mean = run_LDA(X, y)

% Convert labels to categorical
y = categorical(y);

% Ensure classes are stratified
cv = cvpartition(y,'KFold',min(5, min(countcats(y)))); 
% → max K = number of trials in smallest class
% → ensures every fold has at least 1 trial per class

acc = zeros(cv.NumTestSets,1);

for f = 1:cv.NumTestSets
    trainIdx = cv.training(f);
    testIdx  = cv.test(f);

    Xtrain = X(trainIdx,:);
    Xtest  = X(testIdx,:);
    ytrain = y(trainIdx);
    ytest  = y(testIdx);

    % Safety check
    if isempty(Xtrain) || isempty(Xtest)
        warning('Empty fold detected, skipping...');
        continue
    end

    % Fit LDA with regularization
    lda = fitcdiscr(Xtrain, ytrain, 'DiscrimType','linear', ...
                    'Gamma', 0.01, 'Delta',0.01);

    pred = predict(lda, Xtest);
    acc(f) = mean(pred == ytest);
end

acc_mean = mean(acc);

end
