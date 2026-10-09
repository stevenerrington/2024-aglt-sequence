%% Rigorous LDA decoding of identity and ordinal position, from whole-sequence-anchored activity
%
% Companion/contrast to pca_lda_id_position_rigorous.m. That script (like
% the original pca_lda_id_position.m) draws its features from
% sdf_soundAlign_data -- activity re-epoched and re-baselined to EACH
% ELEMENT'S OWN ONSET, then fits a PCA specific to just the 4 stacked
% position-condition epochs. That is a structurally different
% representation from the one behind the manuscript's strongest ordinal
% result: the Fig 2A/B trajectory/linearity finding (frontal PC1 R^2=0.76
% against time-in-sequence) is computed on sdf.sequenceOnset, a single
% CONTINUOUS per-neuron trace aligned to sequence onset and baseline-
% corrected once against the pre-sequence period (see pca_seq_main.m).
% The dominant slow drift that produces that result may simply not be
% recoverable by a PCA fit narrowly on 4 element-onset-locked snapshots.
%
% This script tests position/identity decoding using windows sliced
% directly out of that same continuous, sequence-onset-anchored,
% globally-baselined trace instead -- i.e. it asks whether position is
% decodable from a representation that actually has access to the
% within-trial, cross-element drift dynamics the linearity analysis
% found, rather than one that re-references itself at every element.
%
% Sequence structure (source: stimuli/doc/sequence description.xlsx,
% sheet "Sheet3", cross-checked against pca_seq_main.m's
% (cond_value==seq | cond_value==seq+4) grouping into 4 sequence pairs):
%   cond_value 1,5 -> Grammatical5_novel  : A C G F C
%   cond_value 2,6 -> Grammatical10_novel : A D C G F
%   cond_value 3,7 -> Grammatical_3       : A C F C G
%   cond_value 4,8 -> Grammatical_8       : A D C F C
% Element onset times relative to sequence onset (source:
% stimuli/stimuli_navigator.m, matches Figure 1B): [0 563 1126 1688 2252] ms
%
% Assumes the following are already in the workspace, exactly as set up by
% aglt_analysis_main.m / pca_seq_main.m: spike_log, dirs, auditory_neuron_idx,
% frontal_neuron_idx. Loads data/spike/<session>_<unitDSP>.mat and
% data/processed/<session>.mat (event_table) per neuron, as pca_seq_main.m
% does -- expect this extraction step to take a while over ~1400+ neurons.

clear identity_trials position_trials rigorous_results_wholeseq

rng(20,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1688 2252];  % fixed, rhythmic element onsets (ms) relative to sequence onset
time_win         = 0:413;                   % post-onset analysis window (ms), as in the manuscript
n_time_bins      = 4;                       % non-overlapping bins spanning the 0-413 ms window
n_pcs            = 3;
min_trials       = 3;
n_pseudo_trial   = 15;
k_folds          = 5;
cv_repeats       = 4;
n_boot           = 300;
n_perm           = 300;

bin_edges = round(linspace(0, length(time_win), n_time_bins+1));

seq_letters = { ...
    'A','C','G','F','C'; ...   % cond_value 1 & 5 -- Grammatical5_novel
    'A','D','C','G','F'; ...   % cond_value 2 & 6 -- Grammatical10_novel
    'A','C','F','C','G'; ...   % cond_value 3 & 7 -- Grammatical_3
    'A','D','C','F','C'};      % cond_value 4 & 8 -- Grammatical_8
seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};

identity_classes = {'C','F','G'};                    % @ position 5, as in the manuscript design
position_classes = {'position_2','position_3','position_4','position_5'};  % identity C only

%% ---- Step 1: single-trial, whole-sequence-anchored, baseline z-scored bin features ----
% One struct array entry per neuron; each field is an Ntrials x n_time_bins
% matrix, sliced directly from the continuous sequenceOnset-aligned trace
% at that class's fixed absolute element-onset offset, using ONE global
% pre-sequence baseline for the whole trial (no per-element re-referencing).

n_neurons = size(spike_log,1);
identity_trials(n_neurons,1) = struct();
position_trials(n_neurons,1) = struct();

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
        continue   % missing/unreadable file -- leave this neuron's pools empty, excluded via min_trials
    end

    event_table = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);

    baseline_window = 800:1000;   % -200 to 0 ms pre-sequence, matches pca_seq_main.m
    all_sdf = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all = event_table.cond_value(nonviol_mask);

    baseline_mu_fr  = nanmean(nanmean(all_sdf(:,baseline_window)));
    baseline_std_fr = nanstd(nanmean(all_sdf(:,baseline_window)));

    if isnan(baseline_mu_fr) || isnan(baseline_std_fr) || baseline_std_fr == 0
        continue
    end

    for class_i = 1:length(identity_classes)
        cond_values = lookup_cond_values(seq_letters, seq_cond_values, 5, identity_classes{class_i});
        identity_trials(neuron_i).(identity_classes{class_i}) = extract_bin_features( ...
            all_sdf, cond_value_all, cond_values, element_onset_ms(5), time_win, bin_edges, n_time_bins, ...
            baseline_mu_fr, baseline_std_fr);
    end

    for class_i = 1:length(position_classes)
        pos = class_i + 1;   % position_classes{1} = position_2, i.e. position index 2
        cond_values = lookup_cond_values(seq_letters, seq_cond_values, pos, 'C');
        position_trials(neuron_i).(position_classes{class_i}) = extract_bin_features( ...
            all_sdf, cond_value_all, cond_values, element_onset_ms(pos), time_win, bin_edges, n_time_bins, ...
            baseline_mu_fr, baseline_std_fr);
    end
end

%% ---- Step 2: run the nested, leakage-free decoding for each area x decoding problem ----

decoding_problems = {
    'Identity_Auditory', identity_trials, identity_classes, auditory_neuron_idx;
    'Identity_Frontal',  identity_trials, identity_classes, frontal_neuron_idx;
    'Position_Auditory', position_trials, position_classes, auditory_neuron_idx;
    'Position_Frontal',  position_trials, position_classes, frontal_neuron_idx;
    };

rigorous_results_wholeseq = struct();

for prob_i = 1:size(decoding_problems,1)

    label        = decoding_problems{prob_i,1};
    trial_struct = decoding_problems{prob_i,2};
    classes      = decoding_problems{prob_i,3};
    neuron_idx   = decoding_problems{prob_i,4};

    fprintf('\n=== %s (whole-sequence anchored) ===\n', strrep(label,'_',' '));

    [boot_acc, perm_acc, n_used] = rigorous_lda_decode( ...
        trial_struct, classes, neuron_idx, min_trials, n_pseudo_trial, n_time_bins, ...
        n_pcs, k_folds, cv_repeats, n_boot, n_perm);

    rigorous_results_wholeseq.(label).boot_acc     = boot_acc;
    rigorous_results_wholeseq.(label).perm_acc     = perm_acc;
    rigorous_results_wholeseq.(label).neurons_used = n_used;
    rigorous_results_wholeseq.(label).n_classes    = length(classes);

    obs_acc = median(boot_acc,'omitnan');
    p_perm  = (1 + sum(perm_acc >= obs_acc)) / (1 + sum(~isnan(perm_acc)));

    fprintf('Neurons usable (>= %d real trials/class): %d / %d\n', min_trials, n_used, length(neuron_idx));
    fprintf('Nested-CV accuracy: median = %.3f, 95%% CI [%.3f, %.3f]\n', ...
        obs_acc, prctile(boot_acc,2.5), prctile(boot_acc,97.5));
    fprintf('Theoretical chance = %.3f | Permutation null: median = %.3f, 95th pctile = %.3f\n', ...
        1/length(classes), median(perm_acc,'omitnan'), prctile(perm_acc,95));
    fprintf('Permutation p-value (observed vs null): p = %.4f\n', p_perm);

end

%% ---- Step 3: plot ----

problem_names = fieldnames(rigorous_results_wholeseq);
figure('Renderer','painters','Position',[100 100 700 450]); hold on

bar_x        = 1:length(problem_names);
bar_median   = nan(size(bar_x));
err_lo       = nan(size(bar_x));
err_hi       = nan(size(bar_x));
chance_level = nan(size(bar_x));

for i = 1:length(problem_names)
    acc = rigorous_results_wholeseq.(problem_names{i}).boot_acc;
    bar_median(i)   = median(acc,'omitnan');
    err_lo(i)       = bar_median(i) - prctile(acc,2.5);
    err_hi(i)       = prctile(acc,97.5) - bar_median(i);
    chance_level(i) = 1 / rigorous_results_wholeseq.(problem_names{i}).n_classes;
end

bar(bar_x, bar_median, 0.5, 'FaceColor',[0.4 0.6 0.8]);
errorbar(bar_x, bar_median, err_lo, err_hi, 'k', 'LineStyle','none','LineWidth',1.2)
plot(bar_x, chance_level, 'rd','MarkerFaceColor','r')
set(gca,'XTick',bar_x,'XTickLabel',strrep(problem_names,'_',' '),'XTickLabelRotation',20)
ylabel('Nested cross-validated accuracy')
ylim([0 1]); box off
title('Whole-sequence-anchored LDA decoding (leakage-free, permutation-tested)')
legend({'Median accuracy','95% bootstrap CI','Theoretical chance'},'Location','best')

%% ---- Step 4: area / feature comparisons ----

fprintf('\nAuditory: identity vs position decoding\n');
bootstrap_compare(rigorous_results_wholeseq.Identity_Auditory.boot_acc, rigorous_results_wholeseq.Position_Auditory.boot_acc);

fprintf('\nFrontal: identity vs position decoding\n');
bootstrap_compare(rigorous_results_wholeseq.Identity_Frontal.boot_acc, rigorous_results_wholeseq.Position_Frontal.boot_acc);

fprintf('\nPosition decoding: frontal vs auditory\n');
bootstrap_compare(rigorous_results_wholeseq.Position_Frontal.boot_acc, rigorous_results_wholeseq.Position_Auditory.boot_acc);

fprintf('\nIdentity decoding: frontal vs auditory\n');
bootstrap_compare(rigorous_results_wholeseq.Identity_Frontal.boot_acc, rigorous_results_wholeseq.Identity_Auditory.boot_acc);


%% ================= Local functions =================

function cond_values = lookup_cond_values(seq_letters, seq_cond_values, position, letter)
% Which nonviolation cond_value codes deliver `letter` at ordinal `position`,
% pooling across whichever of the 4 grammatical sequences put it there.
cond_values = [];
for g = 1:size(seq_letters,1)
    if strcmp(seq_letters{g,position}, letter)
        cond_values = [cond_values, seq_cond_values{g}]; %#ok<AGROW>
    end
end
end


function trial_bins = extract_bin_features(all_sdf, cond_value_all, cond_values, onset_ms, time_win, ...
    bin_edges, n_time_bins, baseline_mu_fr, baseline_std_fr)
% Slices the fixed [onset_ms, onset_ms+413] window directly out of the
% continuous sequence-onset-aligned trace for the requested cond_values,
% bins it, and z-scores with the trial's single whole-sequence baseline
% (no per-element re-referencing).

trial_idx = ismember(cond_value_all, cond_values);
if ~any(trial_idx)
    trial_bins = [];
    return
end

cols = 1000 + onset_ms + time_win;   % 1000 = column offset for t=0 (see baseline_window = 800:1000 <-> -200:0 ms)
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


function [boot_acc, perm_acc, n_used] = rigorous_lda_decode(trial_struct, classes, neuron_idx, ...
    min_trials, n_pseudo_trial, n_time_bins, n_pcs, k_folds, cv_repeats, n_boot, n_perm)

n_classes = length(classes);

usable     = true(length(neuron_idx),1);
trial_pool = cell(length(neuron_idx), n_classes);

for ni = 1:length(neuron_idx)
    for ci = 1:n_classes
        v = trial_struct(neuron_idx(ni)).(classes{ci});
        trial_pool{ni,ci} = v;
        if isempty(v) || size(v,1) < min_trials
            usable(ni) = false;
        end
    end
end

use_idx    = find(usable);
n_used     = length(use_idx);
trial_pool = trial_pool(use_idx,:);

boot_acc = nan(n_boot,1);
for boot_i = 1:n_boot
    [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins);
    boot_acc(boot_i) = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats);
end

perm_acc = nan(n_perm,1);
for perm_i = 1:n_perm
    [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins);
    y = y(randperm(length(y)));
    perm_acc(perm_i) = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats);
end

end


function [X, y] = build_pseudo_population(trial_pool, n_pseudo_trial, n_time_bins)

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


function acc = nested_cv_accuracy(X, y, n_pcs, k_folds, cv_repeats)

fold_acc = nan(k_folds*cv_repeats,1);
f = 1;
for rep_i = 1:cv_repeats
    cvp = cvpartition(y, 'KFold', k_folds);
    for fold_i = 1:k_folds
        train_idx = training(cvp, fold_i);
        test_idx  = test(cvp, fold_i);

        mu = mean(X(train_idx,:),1);
        [coeff, scores_train] = pca(X(train_idx,:));
        n_keep = min(n_pcs, size(coeff,2));
        scores_train = scores_train(:,1:n_keep);
        scores_test  = (X(test_idx,:) - mu) * coeff(:,1:n_keep);

        lda_model = fitcdiscr(scores_train, y(train_idx));
        preds = predict(lda_model, scores_test);
        fold_acc(f) = mean(preds == y(test_idx));
        f = f + 1;
    end
end
acc = mean(fold_acc,'omitnan');

end
