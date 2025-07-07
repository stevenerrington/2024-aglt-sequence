clear a_sdf c_sdf d_sdf f_sdf g_sdf pos1_sdf pos2_sdf pos3_sdf pos4_sdf pos5_sdf

for neuron_i = 1:size(spike_log,1)
    valid_trials = cell2mat(sdf_soundAlign_data{neuron_i}(:,9)) & strcmp(sdf_soundAlign_data{neuron_i}(:,5),'nonviol');

    baseline_mu_fr = nanmean(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));
    baseline_std_fr = nanstd(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(strcmp(sdf_soundAlign_data{neuron_i}(:,3),'Baseline'),1))));

    a_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    c_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    d_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    f_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    g_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'),1))),100)-baseline_mu_fr)./baseline_std_fr;

    pos1_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos2_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos3_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos4_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4'),1))),100)-baseline_mu_fr)./baseline_std_fr;
    pos5_sdf(neuron_i,:) = (smooth(nanmean(cell2mat(sdf_soundAlign_data{neuron_i}(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C') & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5'),1))),100)-baseline_mu_fr)./baseline_std_fr;



    n_trl_id(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'A'));
    n_trl_id(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_id(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'D'));
    n_trl_id(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'F'));
    n_trl_id(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'G'));

    n_trl_pos(neuron_i,1) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_1') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,2) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_2') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,3) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_3') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,4) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_4') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));
    n_trl_pos(neuron_i,5) = sum(valid_trials & strcmp(sdf_soundAlign_data{neuron_i}(:,4),'position_5') & strcmp(sdf_soundAlign_data{neuron_i}(:,3),'C'));   

end

time_win = [0:413];

signal_in_identity_auditory = {};  % Initialize a cell array to store the signal data
signal_in_identity_auditory = {c_sdf(auditory_neuron_idx, 200+time_win),...  % SDF for sequence 1 in the current area
    f_sdf(auditory_neuron_idx, 200+time_win),...
    g_sdf(auditory_neuron_idx, 200+time_win)};


signal_in_identity_frontal = {};  % Initialize a cell array to store the signal data
signal_in_identity_frontal = {c_sdf(frontal_neuron_idx, 200+time_win),...  % SDF for sequence 1 in the current area
    f_sdf(frontal_neuron_idx, 200+time_win),...
    g_sdf(frontal_neuron_idx, 200+time_win)};


signal_in_order_auditory = {};  % Initialize a cell array to store the signal data
signal_in_order_auditory = {pos2_sdf(auditory_neuron_idx, 200+time_win),...  % SDF for sequence 1 in the current area
    pos3_sdf(auditory_neuron_idx, 200+time_win),...
    pos4_sdf(auditory_neuron_idx, 200+time_win),...
    pos5_sdf(auditory_neuron_idx, 200+time_win)};


signal_in_order_frontal = {};  % Initialize a cell array to store the signal data
signal_in_order_frontal = {pos2_sdf(frontal_neuron_idx, 200+time_win),...  % SDF for sequence 1 in the current area
    pos3_sdf(frontal_neuron_idx, 200+time_win),...
    pos4_sdf(frontal_neuron_idx, 200+time_win),...
    pos5_sdf(frontal_neuron_idx, 200+time_win)};





% Perform cross-condition PCA analysis to extract principal components
[pc_out_identity_auditory, ~] = get_xcond_pca(signal_in_identity_auditory);  % Perform PCA on the input signals
[pc_out_identity_frontal, ~] = get_xcond_pca(signal_in_identity_frontal);  % Perform PCA on the input signals
[pc_out_order_auditory, ~] = get_xcond_pca(signal_in_order_auditory);  % Perform PCA on the input signals
[pc_out_order_frontal, ~] = get_xcond_pca(signal_in_order_frontal);  % Perform PCA on the input signals




figuren('Renderer', 'painters', 'Position', [200 200 1600 400]);
subplot(1,4,1); hold on
plot3(pc_out_identity_auditory{1}(1,:),pc_out_identity_auditory{1}(2,:),pc_out_identity_auditory{1}(3,:))
plot3(pc_out_identity_auditory{2}(1,:),pc_out_identity_auditory{2}(2,:),pc_out_identity_auditory{2}(3,:))
plot3(pc_out_identity_auditory{3}(1,:),pc_out_identity_auditory{3}(2,:),pc_out_identity_auditory{3}(3,:))
title('ID - Auditory'); view(-45.1291,24.9773); grid on

subplot(1,4,2); hold on
plot3(pc_out_identity_frontal{1}(1,:),pc_out_identity_frontal{1}(2,:),pc_out_identity_frontal{1}(3,:))
plot3(pc_out_identity_frontal{2}(1,:),pc_out_identity_frontal{2}(2,:),pc_out_identity_frontal{2}(3,:))
plot3(pc_out_identity_frontal{3}(1,:),pc_out_identity_frontal{3}(2,:),pc_out_identity_frontal{3}(3,:))
title('ID - Frontal'); view(-45.1291,24.9773); grid on

subplot(1,4,3); hold on
plot3(pc_out_order_auditory{1}(1,:),pc_out_order_auditory{1}(2,:),pc_out_order_auditory{1}(3,:))
plot3(pc_out_order_auditory{2}(1,:),pc_out_order_auditory{2}(2,:),pc_out_order_auditory{2}(3,:))
plot3(pc_out_order_auditory{3}(1,:),pc_out_order_auditory{3}(2,:),pc_out_order_auditory{3}(3,:))
plot3(pc_out_order_auditory{4}(1,:),pc_out_order_auditory{4}(2,:),pc_out_order_auditory{4}(3,:))
title('Position - Auditory'); view(-45.1291,24.9773); grid on

subplot(1,4,4); hold on
plot3(pc_out_order_frontal{1}(1,:),pc_out_order_frontal{1}(2,:),pc_out_order_frontal{1}(3,:))
plot3(pc_out_order_frontal{2}(1,:),pc_out_order_frontal{2}(2,:),pc_out_order_frontal{2}(3,:))
plot3(pc_out_order_frontal{3}(1,:),pc_out_order_frontal{3}(2,:),pc_out_order_frontal{3}(3,:))
plot3(pc_out_order_frontal{4}(1,:),pc_out_order_frontal{4}(2,:),pc_out_order_frontal{4}(3,:))
title('Position - Frontal'); view(-45.1291,24.9773); grid on


%%
clear pc_in
pc_in = pc_out_identity_auditory;

nconds = size(pc_in,2); npcs = 3;

scores = []; labels = [];
for cond_i = 1:nconds
        scores_a = [];
        for pc_i = 1:npcs
            scores_a = [scores_a, pc_in{cond_i}(pc_i,:)'];
        end
        scores = [scores; scores_a];
        labels = [labels, repmat(cond_i,1,length(scores_a))];
end


clear lda_*
lda_model = fitcdiscr(scores, labels);
lda_scores = lda_model.X * lda_model.Coeffs(1,2).Linear;  % For 2-class case

cvmodel = crossval(lda_model, 'KFold', 5);  % 5-fold CV
cv_loss = kfoldLoss(cvmodel);  % classification error
cv_accuracy = 1 - cv_loss;

predicted_labels = predict(lda_model, scores(:,1:npcs));
cm = confusionmat(labels, predicted_labels);
cm_normalized = cm ./ sum(cm, 2);  % Row-wise normalization




figure('Renderer', 'painters', 'Position', [100 100 800 300]); hold on;
subplot(1,2,1)
predicted_labels = predict(lda_model, scores(:,1:npcs));
confusionchart(labels, predicted_labels, 'Normalization', 'row-normalized');
cm = confusionmat(labels, predicted_labels);
cm_normalized = cm ./ sum(cm, 2);  % Row-wise normalization
title('Confusion Matrix');

subplot(1,2,2)
class_accuracy = diag(cm_normalized);
plot(class_accuracy, '-o', 'LineWidth', 2);
xlabel('Class Index');
ylabel('Accuracy');
title('Per-Class Accuracy (Diagonal of Confusion Matrix)');
xlim([0 nconds+1]);
ylim([0 1]);
grid on;
xticks(1:numel(class_accuracy));
xticklabels({'A','C','D','F','G'});  % replace with your actual class labels
hline(0.2,'k--')




























nboot = 1000;
nbootsamples = 25;

area = {'frontal'};

clear pc_out

boot_neurons = [];
for boot_i = 1:nboot
    boot_neurons(:,boot_i) = sort(frontal_neuron_idx(randperm(length(frontal_neuron_idx), nbootsamples)));
end

for boot_i = 1:nboot

    signal_in = {};  % Initialize a cell array to store the signal data
    signal_in = {a_sdf(boot_neurons(:,boot_i), 200+[0:413]),...  % SDF for sequence 1 in the current area
        c_sdf(boot_neurons(:,boot_i), 200+[0:413]),...
        d_sdf(boot_neurons(:,boot_i), 200+[0:413]),...
        f_sdf(boot_neurons(:,boot_i), 200+[0:413]),...
        g_sdf(boot_neurons(:,boot_i), 200+[0:413])};

    % Perform cross-condition PCA analysis to extract principal components
    [pc_out{boot_i}, pc_shuf_out] = get_xcond_pca(signal_in);  % Perform PCA on the input signals


end

%
% for boot_i = 1:10
%     figuren;
%     plot3(pc_out{boot_i}{1}(1,:),pc_out{boot_i}{1}(2,:),pc_out{boot_i}{1}(3,:),'LineWidth',2)
%     plot3(pc_out{boot_i}{2}(1,:),pc_out{boot_i}{2}(2,:),pc_out{boot_i}{2}(3,:),'LineWidth',2)
%     plot3(pc_out{boot_i}{3}(1,:),pc_out{boot_i}{3}(2,:),pc_out{boot_i}{3}(3,:),'LineWidth',2)
%     plot3(pc_out{boot_i}{4}(1,:),pc_out{boot_i}{4}(2,:),pc_out{boot_i}{4}(3,:),'LineWidth',2)
%     plot3(pc_out{boot_i}{5}(1,:),pc_out{boot_i}{5}(2,:),pc_out{boot_i}{5}(3,:),'LineWidth',2)
%     title(int2str(boot_i))
% end

scores = []; labels = [];
npcs = 5;
nconds = 5;

for boot_i = 1:nboot
    scores = []; labels = [];
    for cond_i = 1:nconds
        scores_a = [];
        for pc_i = 1:npcs
            scores_a = [scores_a, pc_out{boot_i}{cond_i}(pc_i,:)'];
        end
        scores = [scores; scores_a];
        labels = [labels, repmat(cond_i,1,length(scores_a))];
    end


    clear lda_*
    lda_model = fitcdiscr(scores, labels);
    lda_scores = lda_model.X * lda_model.Coeffs(1,2).Linear;  % For 2-class case


    cvmodel = crossval(lda_model, 'KFold', 5);  % 5-fold CV
    cv_loss = kfoldLoss(cvmodel);  % classification error
    cv_accuracy = 1 - cv_loss;

    acc_lda_out(boot_i,1) = cv_accuracy;


    predicted_labels = predict(lda_model, scores(:,1:npcs));
    cm = confusionmat(labels, predicted_labels);
    cm_normalized = cm ./ sum(cm, 2);  % Row-wise normalization
    acc_lda_out_ind(boot_i,:) = [diag(cm_normalized)]';
end

labels = [repmat({'A'},nboot,1); repmat({'C'},nboot,1); repmat({'D'},nboot,1); repmat({'F'},nboot,1); repmat({'G'},nboot,1)];

clear lda_accuracy predicted_labels
lda_accuracy(1,1)=gramm('x',labels,'y',[acc_lda_out_ind(:,1); acc_lda_out_ind(:,2); acc_lda_out_ind(:,3); acc_lda_out_ind(:,4); acc_lda_out_ind(:,5)]);
lda_accuracy(1,1).stat_summary('geom',{'line','black_errorbar'},'width',0.5);
% lda_accuracy(1,1).geom_jitter('alpha',0.05);
fig = figure('Renderer', 'painters', 'Position', [100 100 500 400]);
lda_accuracy.draw();



fprintf('Cross-Validated Accuracy: %.2f%%\n', cv_accuracy * 100);

%% 


% Assuming two LDA dimensionsr
Y = lda_model.X * lda_model.Coeffs(1,2).Linear;  % LD1
Z = lda_model.X * lda_model.Coeffs(2,3).Linear;  % LD2


figure;
subplot(1,2,1)
scatter3(scores(:,1), scores(:,2), scores(:,3), 20, categorical(labels), 'filled');
xlabel('LD1'); ylabel('LD2');
title('LDA Projection');

subplot(1,2,2)
scatter3(scores(:,1), scores(:,2), scores(:,3), 20, predicted_labels, 'filled');
xlabel('PC1'); ylabel('PC2');
title('PC Projection');
