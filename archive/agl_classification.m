

for neuron_i = 1:size(spike_log,1)

    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1)); % Display progress for current neuron

    clear a_trials c_trials d_trials f_trials g_trials
    a_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A');
    c_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C');
    d_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D');
    f_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F');
    g_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G');

    % Normalize
    clear baseline_trials bl_fr_mu bl_fr_std
    baseline_trials = strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline');
    bl_fr_mu = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(baseline_trials,1))));
    bl_fr_std = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(baseline_trials,1))));

    a_sdf(neuron_i,:) = (nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(a_trials,1)))-bl_fr_mu)./bl_fr_std;
    c_sdf(neuron_i,:) = (nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(c_trials,1)))-bl_fr_mu)./bl_fr_std;
    d_sdf(neuron_i,:) = (nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(d_trials,1)))-bl_fr_mu)./bl_fr_std;
    f_sdf(neuron_i,:) = (nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(f_trials,1)))-bl_fr_mu)./bl_fr_std;
    g_sdf(neuron_i,:) = (nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(g_trials,1)))-bl_fr_mu)./bl_fr_std;

end

%%

area = 'frontal';

a_class_data = a_sdf(neuron_class.(area).glm_pos,:);
c_class_data = c_sdf(neuron_class.(area).glm_pos,:);
d_class_data = d_sdf(neuron_class.(area).glm_pos,:);
f_class_data = f_sdf(neuron_class.(area).glm_pos,:);
g_class_data = g_sdf(neuron_class.(area).glm_pos,:);

figuren;
plot(ops.sound_sdf_window, nanmean(a_class_data))
plot(ops.sound_sdf_window, nanmean(c_class_data))
plot(ops.sound_sdf_window, nanmean(d_class_data))
plot(ops.sound_sdf_window, nanmean(f_class_data))
plot(ops.sound_sdf_window, nanmean(g_class_data))







%% Accuracy x time

obs_data = {a_class_data, c_class_data, d_class_data, f_class_data, g_class_data};
obs_labels = {'A', 'C', 'D', 'F', 'G'};

accuracy = run_multicat_categorization(obs_data, obs_labels);





%% Figure
figuren;
subplot(5,1,[1 2 3 4]); hold on
plot(window_time, nanmean(window_sdf(strcmp(string(data_labels),'A'),:)))
plot(window_time, nanmean(window_sdf(strcmp(string(data_labels),'C'),:)))
plot(window_time, nanmean(window_sdf(strcmp(string(data_labels),'D'),:)))
plot(window_time, nanmean(window_sdf(strcmp(string(data_labels),'F'),:)))
plot(window_time, nanmean(window_sdf(strcmp(string(data_labels),'G'),:)))


subplot(5,1,5)
plot(accuracy(2,:),accuracy(1,:))
hline(20,'r--')
box off




%% Machine learning core


% Split data into training and testing sets (e.g., 70% training, 30% testing)
cv = cvpartition(data_table.data_labels, 'HoldOut', 0.3);  % 30% test, 70% train
XTrain = data_table.mean_fr(training(cv), :);
YTrain = data_table.data_labels(training(cv), :);
XTest = data_table.mean_fr(test(cv), :);
YTest = data_table.data_labels(test(cv), :);

%Train a Decision Tree classifier
TreeModel = fitctree(XTrain, YTrain);

% Make predictions on the test set
YPred = predict(TreeModel, XTest);

% Calculate accuracy
accuracy = sum(YPred == YTest) / length(YTest);
fprintf('Accuracy: %.2f%%\n', accuracy * 100);

% Display confusion matrix
figure;
confusionchart(YTest, YPred);


