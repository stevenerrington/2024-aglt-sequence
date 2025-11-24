%% ===============================
% TRIAL-BASED PCA + LDA
% ===============================

time_window = [0 563];   % ms window
window_idx = find(ops.sound_sdf_window >= time_window(1) & ops.sound_sdf_window <= time_window(2));

X_id_raw = [];  % trial x neuron
y_id     = [];  % trial labels

X_pos_raw = [];
y_pos     = [];

% -------------------------------
% Step 1: Loop over neurons and extract trial-mean firing rates
% -------------------------------
for neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i/%i\n', neuron_i, size(spike_log,1));
    trial_sdf     = sdf_soundAlign_data{neuron_i}(:,1);
    trial_id      = sdf_soundAlign_data{neuron_i}(:,3);
    trial_pos     = sdf_soundAlign_data{neuron_i}(:,4);
    trial_type    = sdf_soundAlign_data{neuron_i}(:,5);
    trial_correct = sdf_soundAlign_data{neuron_i}(:,9);

    valid_mask = cell2mat(trial_correct) & strcmp(trial_type,'nonviol');

    % ---- Identity (position 5 only) ----
    id_mask = valid_mask & strcmp(trial_pos,'position_5');
    for t = find(id_mask')
        X_id_raw(end+1, 1) = mean(trial_sdf{t}(window_idx));
        y_id{end+1,1} = trial_id{t};
        y2_id{end+1,1} = spike_log.area{neuron_i};
    end

    % ---- Position (identity C only) ----
    pos_mask = valid_mask & strcmp(trial_id,'C');
    for t = find(pos_mask')
        X_pos_raw(end+1, 1) = mean(trial_sdf{t}(window_idx));
        y_pos{end+1,1} = trial_pos{t};
        y2_pos{end+1,1} = spike_log.area{neuron_i};
    end
end

y_id  = categorical(y_id);
y_pos = categorical(y_pos);

% -------------------------------
% Step 2: PCA for dimensionality reduction
% -------------------------------
nPCs = 5;  % number of PCs to keep

% Identity
[coeff_id, score_id, ~, ~, explained_id] = pca(X_id_raw);
X_id = score_id(:,1:nPCs);

% Position
[coeff_pos, score_pos, ~, ~, explained_pos] = pca(X_pos_raw);
X_pos = score_pos(:,1:nPCs);

% -------------------------------
% Step 3: Run LDA on trial-wise PCA scores
% -------------------------------
acc_id  = run_LDA(X_id, y_id);
acc_pos = run_LDA(X_pos, y_pos);