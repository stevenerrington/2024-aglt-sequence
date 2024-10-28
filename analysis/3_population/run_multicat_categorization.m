function accuracy = run_multicat_categorization(obs_data, obs_labels)

%% Restructure data
n_cats = length(obs_labels);
data_in = []; labels_in = [];

for cat_i = 1:n_cats
    data_in = [data_in; obs_data{cat_i}];
    labels_in = [labels_in; repmat(obs_labels(cat_i),size(obs_data{cat_i},1),1)];
end

labels_in = categorical(labels_in);

%% Smooth data (for SDFs)
window_size = 50;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(data_in, window_size, window_shift);

%% Generate a data table to pull from
data_table = table(window_sdf, labels_in);

%% Setup cross-validation parameters
cv = cvpartition(data_table.labels_in, 'HoldOut', 0.2);  % 20% test, 80% train
training_idx = training(cv);
test_idx = test(cv);

accuracy = [];

for time_i = 1:length(window_time)
    clear XTrain YTrain XTest YTest TreeModel YPred

    % Split data into training and testing sets (e.g., 70% training, 30% testing)
    XTrain = data_table.window_sdf(training_idx, time_i);
    YTrain = data_table.labels_in(training_idx, :);
    XTest = data_table.window_sdf(test_idx, time_i);
    YTest = data_table.labels_in(test_idx, :);

    %Train a Decision Tree classifier
    TreeModel = fitctree(XTrain, YTrain);

    % Make predictions on the test set
    YPred = predict(TreeModel, XTest);

    % Calculate accuracy
    accuracy(1,time_i) = (sum(YPred == YTest) / length(YTest))*100;

end

accuracy(2,:) = window_time;
