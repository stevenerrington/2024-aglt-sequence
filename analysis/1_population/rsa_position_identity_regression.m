%% RSA multiple regression: does population activity encode POSITION or
%  IDENTITY, once SEQUENCE-GROUP is partialled out in the same model?
%
% Everything so far (position_crossgroup_generalization.m,
% pca_seq_linearity_crossgroup.m) tested one specific comparison at a
% time using classification/regression accuracy, and had to work around
% the position/identity <-> sequence-group confound by picking a single
% train/test split. This script is a genuinely different approach, built
% from scratch rather than as a patch to the classification pipeline:
% representational similarity analysis (RSA) with multiple regression,
% which asks how much UNIQUE variance in population activity patterns is
% explained by "same position", "same identity" and "same sequence-group"
% SIMULTANEOUSLY, rather than testing one factor at a time.
%
% Design. Excluding position 1 (always sound "A" in every sequence, so
% uninformative), the 4 sequences x 4 positions (2-5) give 16 distinct
% (group, position) cells, each tied to a single real physical sequence
% -- no pooling across sequences within a cell, unlike earlier scripts.
% These 16 cells cover all 11 unique (letter, position) combinations that
% exist in this stimulus set (some letters recur at a given position via
% a different sequence, e.g. "C" occurs at position 2 via sequence 1 AND
% via sequence 3 -- both are included as separate cells).
%
% Per bootstrap iteration: draw a pseudo-population pseudo-trial for each
% of the 16 cells (independent per-neuron resampling from each neuron's
% own real trial pool for that cell, as in every earlier script), then
% compute the full pairwise correlation-distance matrix between all
% pseudo-trials (1 - Pearson r across the neuron x time-bin feature
% vector -- standard RSA pattern-similarity). Three same/different
% predictors are built for every pseudo-trial pair: same_position,
% same_identity, same_group. A single OLS regression,
%   distance ~ 1 + same_position + same_identity + same_group
% is fit on ALL pairs at once, so the coefficient on same_position is the
% variance in pattern similarity explained by matching position AFTER
% controlling for matching identity and matching sequence-group (and
% vice versa for the other two coefficients) -- exactly the "partial out
% the confound in the same model" approach described as the fresh
% alternative to the classification pipeline.
%
% Sign convention: population activity that genuinely encodes a variable
% should make same-label pairs MORE similar (smaller distance) than
% different-label pairs, so a real effect shows up as a NEGATIVE
% coefficient. This is called out again in the console output.
%
% Significance: pairwise distances aren't independent observations (each
% pseudo-trial appears in many pairs), so significance is assessed by
% permuting which (group, position, identity) label-tuple is attached to
% which of the 16 underlying data blocks (not permuting individual pairs,
% which would violate the block structure), keeping the real, unshuffled
% distance matrix fixed, and refitting the same regression -- giving a
% paired null coefficient from the exact same iteration's neural draw,
% in the same style as pca_seq_linearity_crossgroup.m's shuffled-null.
%
% Scope note: this is deliberately the "core" version -- one fixed
% 0-413 ms window per position (matching every earlier script), no
% previous-trial-carryover term, and pooled across sessions rather than
% checked session-by-session. Those are natural follow-on extensions,
% not run here.
%
% Assumes the following are already in the workspace: spike_log, dirs,
% auditory_neuron_idx, frontal_neuron_idx.

clear cell_trials rsa_results

rng(20,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1688 2252];
time_win         = 0:413;
n_time_bins      = 4;
min_trials       = 3;
n_pseudo_trial   = 10;   % pseudo-trials per (group,position) cell, per iteration
n_boot           = 300;

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));

seq_letters = { ...
    'A','C','G','F','C'; ...   % group 1 (cond_value 1,5)
    'A','D','C','G','F'; ...   % group 2 (cond_value 2,6)
    'A','C','F','C','G'; ...   % group 3 (cond_value 3,7)
    'A','D','C','F','C'};      % group 4 (cond_value 4,8)
seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};
letter_alphabet = {'A','C','D','F','G'};

% ---- Build the 16 (group,position) cells, positions 2-5 ----
cell_group    = nan(1,16);
cell_position = nan(1,16);
cell_letter   = nan(1,16);
cell_field    = cell(1,16);

cell_i = 0;
for g = 1:4
    for p = 2:5
        cell_i = cell_i + 1;
        cell_group(cell_i)    = g;
        cell_position(cell_i) = p;
        cell_letter(cell_i)   = find(strcmp(letter_alphabet, seq_letters{g,p}));
        cell_field{cell_i}    = sprintf('g%d_p%d', g, p);
    end
end

fprintf('16 (group,position) cells -> letter:\n');
for cell_i = 1:16
    fprintf('  %s : group %d, position %d, letter %s\n', cell_field{cell_i}, ...
        cell_group(cell_i), cell_position(cell_i), letter_alphabet{cell_letter(cell_i)});
end

%% ---- Step 1: per-neuron, per-cell single-trial binned features ----

n_neurons = size(spike_log,1);
cell_trials(n_neurons,1) = struct();

for neuron_i = 1:n_neurons

    if mod(neuron_i,100) == 0
        fprintf('Neuron %i of %i\n', neuron_i, n_neurons);
    end

    try
        sdf_in = load(fullfile(dirs.root, 'data', 'spike', ...
            [spike_log.session{neuron_i} '_' spike_log.unitDSP{neuron_i} '.mat']));
        event_table_in = load(fullfile(dirs.mat_data, ...
            [spike_log.session{neuron_i} '.mat']), 'event_table');
    catch
        continue
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    baseline_window = 800:1000;
    all_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all = event_table.cond_value(nonviol_mask);

    baseline_mu_fr  = nanmean(nanmean(all_sdf(:,baseline_window)));
    baseline_std_fr = nanstd(nanmean(all_sdf(:,baseline_window)));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for cell_i = 1:16
        cond_values = seq_cond_values{cell_group(cell_i)};
        cell_trials(neuron_i).(cell_field{cell_i}) = extract_bin_features( ...
            all_sdf, cond_value_all, cond_values, element_onset_ms(cell_position(cell_i)), ...
            time_win, bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr);
    end
end

%% ---- Step 2: bootstrap RSA regression, per area ----

decoding_problems = {'Auditory', auditory_neuron_idx; 'Frontal', frontal_neuron_idx};

rsa_results = struct();

for prob_i = 1:size(decoding_problems,1)

    label      = decoding_problems{prob_i,1};
    neuron_idx = decoding_problems{prob_i,2};

    usable     = true(length(neuron_idx),1);
    trial_pool = cell(length(neuron_idx), 16);

    for ni = 1:length(neuron_idx)
        for cell_i = 1:16
            v = cell_trials(neuron_idx(ni)).(cell_field{cell_i});
            trial_pool{ni,cell_i} = v;
            if isempty(v) || size(v,1) < min_trials
                usable(ni) = false;
            end
        end
    end

    use_idx    = find(usable);
    n_used     = length(use_idx);
    trial_pool = trial_pool(use_idx,:);

    fprintf('\n=== %s: RSA multiple regression ===\n', label);
    fprintf('Neurons usable (>=%d real trials in all 16 cells): %d / %d\n', min_trials, n_used, length(neuron_idx));

    n_total = 16 * n_pseudo_trial;
    mask    = triu(true(n_total),1);
    n_pairs = sum(mask(:));

    beta_true = nan(n_boot,3);   % columns: position, identity, group
    beta_null = nan(n_boot,3);

    position_full = repelem(cell_position, n_pseudo_trial)';
    letter_full   = repelem(cell_letter,   n_pseudo_trial)';
    group_full    = repelem(cell_group,    n_pseudo_trial)';

    for boot_i = 1:n_boot

        [X, ~] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins);

        R = corrcoef(X');
        D = 1 - R;
        d_vec = D(mask);

        % True labels
        pos_mat = position_full == position_full';
        id_mat  = letter_full   == letter_full';
        grp_mat = group_full    == group_full';
        Xdes = [ones(n_pairs,1), double(pos_mat(mask)), double(id_mat(mask)), double(grp_mat(mask))];
        b = Xdes \ d_vec;
        beta_true(boot_i,:) = b(2:4)';

        % Permuted labels: shuffle which (group,position,letter) tuple
        % belongs to which of the 16 data blocks; same neural distances.
        perm_order = randperm(16);
        position_perm = repelem(cell_position(perm_order), n_pseudo_trial)';
        letter_perm   = repelem(cell_letter(perm_order),   n_pseudo_trial)';
        group_perm    = repelem(cell_group(perm_order),    n_pseudo_trial)';
        pos_mat_p = position_perm == position_perm';
        id_mat_p  = letter_perm   == letter_perm';
        grp_mat_p = group_perm    == group_perm';
        Xdes_p = [ones(n_pairs,1), double(pos_mat_p(mask)), double(id_mat_p(mask)), double(grp_mat_p(mask))];
        b_null = Xdes_p \ d_vec;
        beta_null(boot_i,:) = b_null(2:4)';

    end

    rsa_results.(label).beta_true = beta_true;
    rsa_results.(label).beta_null = beta_null;
    rsa_results.(label).n_used    = n_used;

    predictor_names = {'same_position','same_identity','same_group'};
    fprintf('(negative beta = same-label pairs are MORE similar than different-label pairs, i.e. evidence of encoding)\n');
    for pred_i = 1:3
        fprintf('%s -- true beta: median = %.4f, 95%% CI [%.4f, %.4f] | null beta: median = %.4f\n', ...
            predictor_names{pred_i}, median(beta_true(:,pred_i),'omitnan'), ...
            prctile(beta_true(:,pred_i),2.5), prctile(beta_true(:,pred_i),97.5), ...
            median(beta_null(:,pred_i),'omitnan'));
        fprintf('  %s: true vs null:\n', predictor_names{pred_i});
        bootstrap_compare(beta_true(:,pred_i), beta_null(:,pred_i));
    end

end

%% ---- Step 3: plot ----

figure('Renderer','painters','Position',[100 100 550 450]); hold on

predictor_names = {'same position','same identity','same group'};
areas = {'Auditory','Frontal'};
bar_x = 1:6;
bar_median = nan(1,6); err_lo = nan(1,6); err_hi = nan(1,6);
col_i = 1;
for pred_i = 1:3
    for a = 1:2
        v = rsa_results.(areas{a}).beta_true(:,pred_i);
        bar_median(col_i) = median(v,'omitnan');
        err_lo(col_i) = bar_median(col_i) - prctile(v,2.5);
        err_hi(col_i) = prctile(v,97.5) - bar_median(col_i);
        col_i = col_i + 1;
    end
end

bar(bar_x, bar_median, 0.6, 'FaceColor',[0.4 0.6 0.8]);
errorbar(bar_x, bar_median, err_lo, err_hi, 'k', 'LineStyle','none','LineWidth',1.2)
plot(bar_x, zeros(size(bar_x)), 'k--')
set(gca,'XTick',bar_x,'XTickLabel',{'Aud pos','Frontal pos','Aud identity','Frontal identity','Aud group','Frontal group'},'XTickLabelRotation',30)
ylabel('Regression coefficient (distance ~ same-label)')
box off
title('RSA multiple regression: unique contribution of position/identity/group')

%% ---- Step 4: frontal vs auditory, on the position and identity coefficients ----

fprintf('\nPosition coefficient: frontal vs auditory (more negative = stronger position encoding)\n');
bootstrap_compare(rsa_results.Frontal.beta_true(:,1), rsa_results.Auditory.beta_true(:,1));

fprintf('\nIdentity coefficient: frontal vs auditory\n');
bootstrap_compare(rsa_results.Frontal.beta_true(:,2), rsa_results.Auditory.beta_true(:,2));


%% ================= Local functions =================

function trial_bins = extract_bin_features(all_sdf, cond_value_all, cond_values, onset_ms, time_win, ...
    bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr)
% Identical to earlier scripts' version.

trial_idx = ismember(cond_value_all, cond_values);
if ~any(trial_idx)
    trial_bins = [];
    return
end

cols = 1000 + onset_ms + time_win;
trial_window = all_sdf(trial_idx, cols);

raw_bins = nan(size(trial_window,1), n_time_bins);
for bin_i = 1:n_time_bins
    raw_bins(:,bin_i) = nanmean(trial_window(:, bin_edges(bin_i)+1:bin_edges(bin_i+1)), 2);
end

valid_rows = all(~isnan(raw_bins),2);
raw_bins = raw_bins(valid_rows,:);

if isempty(raw_bins)
    trial_bins = [];
else
    trial_bins = (raw_bins - baseline_mu_fr) ./ baseline_std_fr;
end

end


function [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins)
% Identical to earlier scripts' version -- generic over however many
% classes trial_pool has columns for (16, here).

[n_neurons, n_classes] = size(trial_pool);
X = nan(n_pseudo_trial*n_classes, n_neurons*n_time_bins);
y = nan(n_pseudo_trial*n_classes, 1);

row0 = 0;
for ci = 1:n_classes
    block = nan(n_pseudo_trial, n_neurons*n_time_bins);
    for ni = 1:n_neurons
        pool = trial_pool{ni,ci};
        draw_idx = randi(size(pool,1), n_pseudo_trial, 1);
        cols = (ni-1)*n_time_bins + (1:n_time_bins);
        block(:,cols) = pool(draw_idx,:);
    end
    X(row0+1:row0+n_pseudo_trial,:) = block;
    y(row0+1:row0+n_pseudo_trial)   = ci;
    row0 = row0 + n_pseudo_trial;
end

end
