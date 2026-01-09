function [acc_identity, acc_position, details] = decode_identity_position_singleNeuron(neuron_data)
%--------------------------------------------------------------------------
% Decode sound identity and sequence position for a single neuron
%
% Uses:
%   - valid trials only
%   - nonviol trials only
%   - excludes Baseline and A
%   - balanced classes
%   - linear SVM with cross-validation
%--------------------------------------------------------------------------

%% Parameters

time_vec   = -200:800;
decode_win = time_vec >= 0 & time_vec <= 413;
nFolds     = 5;
nPerm = 1000;

details = struct();

%% Trial filtering

is_valid   = cell2mat(neuron_data(:,9)) == 1;
is_nonviol = strcmp(neuron_data(:,5), 'nonviol');
is_sound   = ~strcmp(neuron_data(:,3), 'Baseline');
not_A      = ~strcmp(neuron_data(:,3), 'A');

keep_idx = is_valid & is_nonviol & is_sound & not_A;
data = neuron_data(keep_idx,:);

if size(data,1) < 10
    acc_identity = NaN;
    acc_position = NaN;
    return
end

%% Feature extraction (mean FR)

nTrials = size(data,1);
X = nan(nTrials,1);

for t = 1:nTrials
    sdf = data{t,1};
    X(t) = mean(sdf(decode_win));
end

%% =========================
% Identity decoding
% =========================

Y_id = categorical(data(:,3));
id_counts = countcats(Y_id);

if numel(id_counts) < 2 || min(id_counts) < nFolds
    acc_identity = NaN;
else
    [X_id, Y_id] = balance_classes(X, Y_id);

    [acc_identity, confMat_id, classAcc_id] = ...
        svm_decode_cv(X_id, Y_id, nFolds);

    details.identity.labels   = categories(Y_id);
    details.identity.classAcc = classAcc_id;
    details.identity.confMat  = confMat_id;
end


[chance_id, null_id] = ...
    svm_decode_permutation(X_id, Y_id, nFolds, nPerm);

details.identity.chanceThr = chance_id;
details.identity.nullAcc   = null_id;


%% =========================
% Position decoding
% =========================

Y_pos = categorical(data(:,4));

% Correct categorical comparison
valid_pos = Y_pos ~= categorical("position_1");

X_pos = X(valid_pos,:);
Y_pos = Y_pos(valid_pos);

pos_counts = countcats(Y_pos);

if numel(pos_counts) < 2 || min(pos_counts) < nFolds
    acc_position = NaN;
else
    [X_pos, Y_pos] = balance_classes(X_pos, Y_pos);

    [acc_position, confMat_pos, classAcc_pos] = ...
        svm_decode_cv(X_pos, Y_pos, nFolds);

    details.position.labels   = categories(Y_pos);
    details.position.classAcc = classAcc_pos;
    details.position.confMat  = confMat_pos;
end


[chance_pos, null_pos] = ...
    svm_decode_permutation(X_pos, Y_pos, nFolds, nPerm);

details.position.chanceThr = chance_pos;
details.position.nullAcc   = null_pos;

end
