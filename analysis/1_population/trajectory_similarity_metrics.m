%% Trajectory-similarity metrics in full neuron space (no PCA), cross-validated
%% against trial noise, with neuron-resampling CIs and smoothness-safe nulls
%
% WHAT ALREADY EXISTS (pca_seq_main.m): path length, curvature, net
% displacement, efficiency, radius of gyration, tortuosity, hull volume, mean
% cosine similarity, element-onset distance, start-end distance, and Procrustes
% similarity between consecutive elements. They are computed in a 3-PC space
% of the trial-AVERAGED SDF, with a neuron bootstrap (500 neurons drawn with
% replacement) and no noise floor or null. Procrustes also removes translation
% and rotation, so it measures within-element SHAPE only and cannot tell
% "state reset" from "state advance".
%
% WHAT THIS SCRIPT ADDS
% Each epoch (0-413 ms after an element onset) is a vector of baseline-z-scored
% firing rates (neurons x time bins). Distances and similarities between
% epochs are CROSS-VALIDATED: epoch means come from two disjoint trial halves
% (or, in Part 2, from two disjoint sequence groups) and the two halves are
% multiplied, so independent noise cancels in expectation. The result is an
% unbiased estimate of the squared distance (or cosine) with no noise floor to
% subtract. Because every statistic is a mean of per-neuron terms, resampling
% NEURONS (the real sampling unit; neurons are not simultaneously recorded) is
% exact and cheap.
%
%   PART 1  Descriptive trajectory geometry (what the paper's trajectories
%           show): all sequences pooled, epochs 1-5.
%       - cosine similarity between epochs (1 = same population state)
%       - relative dissimilarity D(k,j) / mean response energy
%       - "reset index" (mean over all epoch pairs), adjacent-epoch value, and
%         DRIFT WITH LAG: Spearman correlation and slope of dissimilarity vs
%         |k-j|. A resetting code is flat in lag; a progressing code grows.
%       - Null for the lag relation: exact permutation of epoch ORDER (all
%         5! = 120 orders applied to every neuron together), which preserves
%         smoothness within epochs. Frontal-vs-auditory by neuron bootstrap.
%       CAVEAT: epochs differ in element identity as well as position, so this
%       describes how the state moves, not what it encodes.
%
%   PART 2  Identity-matched, sequence-controlled version: letter C at
%           positions 2, 3, 5. Sequence groups: A = {1,2}, B = {3,4}.
%       - CROSS-GROUP distance: (mean_A(p) - mean_A(q)) . (mean_B(p) - mean_B(q)).
%         Only the part of the position difference that is present in BOTH
%         sequence groups survives (sequence-specific parts are independent
%         across groups and cancel in expectation).
%       - WITHIN-GROUP distance: same quantity from random trial halves that
%         both contain both groups (position difference INCLUDING sequence
%         context).
%       - Fraction sequence-independent = cross-group / within-group.
%       - Null: per-neuron shuffling of position labels among trials within
%         each sequence group (500 permutations), recomputed for every neuron
%         and averaged by area.
%
% PRE-SPECIFIED READING (written before any result exists)
%   "Frontal trajectories progress while auditory ones reset" is supported only
%   if, in Part 1, frontal dissimilarity rises with lag (Spearman rho > 0,
%   epoch-order permutation p < 0.05) AND the frontal lag slope exceeds
%   auditory (neuron bootstrap p < 0.05).
%   "A sequence-independent position-related signal, stronger in frontal" is
%   supported only if, in Part 2, the frontal cross-group distance exceeds its
%   permutation null (p < 0.05) AND exceeds auditory (bootstrap and permutation
%   p < 0.05). Part 1 alone is trajectory geometry, not a position code.
%   Position is perfectly confounded with elapsed time in this design, so a
%   positive Part 2 means "position / elapsed time that generalises across
%   sequences".
%
% Assumes in the workspace: spike_log, dirs, auditory_neuron_idx,
% frontal_neuron_idx (as set up by aglt_analysis_main.m).
%
% Alignment notes: element onsets use the main-pipeline values
% [0 563 1126 1689 2252]; the sdf column for time t is 1001 + t (timewin starts
% at -1000), so zero_col = 1001 here.

clear Pn1 DXn SXn DWn SWn DXnull

rng(41,'twister')

%% ---- User-set parameters ----
element_onset_ms = [0 563 1126 1689 2252];
time_win         = 0:413;
n_bins           = 8;
bin_edges        = round(linspace(0, length(time_win), n_bins+1));
zero_col         = 1001;

n_splits         = 20;      % random trial split-halves averaged per neuron
n_boot           = 2000;    % neuron-resampling bootstrap
n_perm           = 500;     % Part 2 label permutations
min_trials_epoch = 20;      % Part 1: trials per neuron
min_trials_c     = 3;       % Part 2: trials per (position x sequence group) cell
baseline_window  = 800:1000;

seq_letters = { ...
    'A','C','G','F','C'; ...   % cond_value 1 & 5   (group A)
    'A','D','C','G','F'; ...   % cond_value 2 & 6   (group A)
    'A','C','F','C','G'; ...   % cond_value 3 & 7   (group B)
    'A','D','C','F','C'};      % cond_value 4 & 8   (group B)
seq_cond_values = {[1 5],[2 6],[3 7],[4 8]};
groupset = {[1 2],[3 4]};
pos_c    = [2 3 5];
pairs    = [1 2; 2 3; 1 3];     % position pairs in Part 2 (indices into pos_c)

cols_area = [130 3 51; 41 107 115]./255;   % auditory, frontal (project palette)
area_names = {'Auditory','Frontal'};

%% ---- Step 1: per-neuron cross-validated terms ----

n_neurons = size(spike_log,1);
Pn1   = nan(n_neurons,5,5);
DXn   = nan(n_neurons,3);
SXn   = nan(n_neurons,3);
DWn   = nan(n_neurons,3);
SWn   = nan(n_neurons,3);
DXnull = nan(n_neurons,n_perm,3,'single');

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

    event_table  = event_table_in.event_table;
    nonviol_mask = strcmp(event_table.cond_label,'nonviol') & ~isnan(event_table.rewardOnset_ms);
    all_sdf        = sdf_in.sdf.sequenceOnset(nonviol_mask,:);
    cond_value_all = event_table.cond_value(nonviol_mask);

    bl_trial = mean(all_sdf(:,baseline_window), 2, 'omitnan');
    mu_fr = mean(bl_trial,'omitnan');
    sd_fr = std(bl_trial,'omitnan');
    if isnan(mu_fr) || isnan(sd_fr) || sd_fr == 0
        continue
    end

    % ---- Part 1: all sequences pooled, five epochs ----
    nt = size(all_sdf,1);
    F  = nan(nt,5,n_bins);
    for e = 1:5
        w = all_sdf(:, zero_col + element_onset_ms(e) + time_win);
        for b = 1:n_bins
            F(:,e,b) = mean(w(:, bin_edges(b)+1:bin_edges(b+1)), 2, 'omitnan');
        end
    end
    F = (F - mu_fr) ./ sd_fr;
    good = all(all(~isnan(F),3),2);
    F = F(good,:,:);
    ng = size(F,1);
    if ng >= min_trials_epoch
        Psum = zeros(5,5);
        for r = 1:n_splits
            idx = randperm(ng);
            h   = floor(ng/2);
            m1  = reshape(mean(F(idx(1:h),:,:),1), 5, n_bins);
            m2  = reshape(mean(F(idx(h+1:2*h),:,:),1), 5, n_bins);
            Psum = Psum + 0.5*(m1*m2' + m2*m1');
        end
        Pn1(neuron_i,:,:) = Psum / n_splits;
    end

    % ---- Part 2: letter C at positions 2,3,5, by sequence group ----
    X = cell(3,2);
    ok = true;
    for pi_ = 1:3
        for g = 1:2
            cv = [];
            for sg = groupset{g}
                if strcmp(seq_letters{sg, pos_c(pi_)}, 'C')
                    cv = [cv, seq_cond_values{sg}]; %#ok<AGROW>
                end
            end
            X{pi_,g} = cell_feats(all_sdf, cond_value_all, cv, element_onset_ms(pos_c(pi_)), ...
                zero_col, time_win, bin_edges, n_bins, mu_fr, sd_fr);
            if isempty(X{pi_,g}) || size(X{pi_,g},1) < min_trials_c
                ok = false;
            end
        end
    end
    if ~ok; continue; end

    % cross-group
    mA = zeros(3,n_bins); mB = zeros(3,n_bins);
    for pi_ = 1:3
        mA(pi_,:) = mean(X{pi_,1},1);
        mB(pi_,:) = mean(X{pi_,2},1);
    end
    PX = 0.5*(mA*mB' + mB*mA');
    SX = diag(PX)';
    SXn(neuron_i,:) = SX;
    for pr = 1:3
        p = pairs(pr,1); q = pairs(pr,2);
        DXn(neuron_i,pr) = SX(p) + SX(q) - 2*PX(p,q);
    end

    % within-group (random halves stratified by cell; sequence context NOT removed)
    PWs = zeros(3,3);
    for r = 1:n_splits
        m1 = zeros(3,n_bins); m2 = zeros(3,n_bins);
        for pi_ = 1:3
            t1 = []; t2 = [];
            for g = 1:2
                v = X{pi_,g};
                idx = randperm(size(v,1));
                h = floor(size(v,1)/2);
                t1 = [t1; v(idx(1:h),:)];             %#ok<AGROW>
                t2 = [t2; v(idx(h+1:2*h),:)];         %#ok<AGROW>
            end
            m1(pi_,:) = mean(t1,1);
            m2(pi_,:) = mean(t2,1);
        end
        PWs = PWs + 0.5*(m1*m2' + m2*m1');
    end
    PW = PWs / n_splits;
    SW = diag(PW)';
    SWn(neuron_i,:) = SW;
    for pr = 1:3
        p = pairs(pr,1); q = pairs(pr,2);
        DWn(neuron_i,pr) = SW(p) + SW(q) - 2*PW(p,q);
    end

    % null: shuffle position labels among trials within each sequence group
    pooled = cell(1,2); lab0 = cell(1,2);
    for g = 1:2
        pooled{g} = [X{1,g}; X{2,g}; X{3,g}];
        lab0{g}   = [ones(size(X{1,g},1),1); 2*ones(size(X{2,g},1),1); 3*ones(size(X{3,g},1),1)];
    end
    for perm_i = 1:n_perm
        mAp = zeros(3,n_bins); mBp = zeros(3,n_bins);
        lA = lab0{1}(randperm(numel(lab0{1})));
        lB = lab0{2}(randperm(numel(lab0{2})));
        for pi_ = 1:3
            mAp(pi_,:) = mean(pooled{1}(lA==pi_,:),1);
            mBp(pi_,:) = mean(pooled{2}(lB==pi_,:),1);
        end
        PXp = 0.5*(mAp*mBp' + mBp*mAp');
        SXp = diag(PXp)';
        for pr = 1:3
            p = pairs(pr,1); q = pairs(pr,2);
            DXnull(neuron_i,perm_i,pr) = SXp(p) + SXp(q) - 2*PXp(p,q);
        end
    end
end

area_idx = {auditory_neuron_idx, frontal_neuron_idx};

%% ---- Step 2: Part 1 (trajectory geometry, all sequences pooled) ----

fprintf('\n================ PART 1: trajectory geometry, all sequences pooled (epochs 1-5) ================\n');
orders = perms(1:5);                 % 120 epoch orders for the lag permutation
boot1 = struct();
obs1  = cell(1,2);
for a = 1:2
    idx = area_idx{a};
    idx = idx(~isnan(Pn1(idx,1,1)));
    N1(a) = numel(idx);
    Pmean = reshape(mean(Pn1(idx,:,:),1), 5, 5);
    obs1{a} = traj_metrics(Pmean);

    % neuron bootstrap
    b_reset = nan(n_boot,1); b_adj = nan(n_boot,1); b_rho = nan(n_boot,1); b_slope = nan(n_boot,1);
    b_cosadj = nan(n_boot,1); b_cosend = nan(n_boot,1); b_lag = nan(n_boot,4);
    for b = 1:n_boot
        bi = idx(randi(numel(idx), numel(idx), 1));
        mm = traj_metrics(reshape(mean(Pn1(bi,:,:),1), 5, 5));
        b_reset(b) = mm.reset; b_adj(b) = mm.adj; b_rho(b) = mm.rho; b_slope(b) = mm.slope;
        b_cosadj(b) = mm.cos_adj; b_cosend(b) = mm.cos_startend; b_lag(b,:) = mm.lagrel;
    end
    boot1(a).reset = b_reset; boot1(a).adj = b_adj; boot1(a).rho = b_rho; boot1(a).slope = b_slope;
    boot1(a).cos_adj = b_cosadj; boot1(a).cos_end = b_cosend; boot1(a).lag = b_lag;

    % exact epoch-order permutation test for dissimilarity-vs-lag
    [pk, pj] = find(triu(true(5),1));
    relv = obs1{a}.rel(sub2ind([5 5], pk, pj));
    rho_perm = nan(size(orders,1),1);
    for oi = 1:size(orders,1)
        pos_of = orders(oi,:);
        lagp = abs(pos_of(pk) - pos_of(pj));
        rho_perm(oi) = corr(lagp(:), relv(:), 'Type','Spearman', 'Rows','complete');
    end
    p_order = mean(abs(rho_perm) >= abs(obs1{a}.rho) - 1e-12);

    fprintf('\n%s (%d neurons)\n', area_names{a}, N1(a));
    fprintf('  cosine similarity between epochs (cross-validated; 1 = same state):\n');
    disp(round(obs1{a}.cs, 2));
    fprintf('  relative dissimilarity D/energy (0 = same state):\n');
    disp(round(obs1{a}.rel, 2));
    fprintf('  reset index (mean dissimilarity over all epoch pairs) = %.3f  [%.3f, %.3f]\n', ...
        obs1{a}.reset, prctile(b_reset,2.5), prctile(b_reset,97.5));
    fprintf('  adjacent-epoch dissimilarity = %.3f  [%.3f, %.3f]; adjacent cosine = %.3f  [%.3f, %.3f]; epoch1-5 cosine = %.3f  [%.3f, %.3f]\n', ...
        obs1{a}.adj, prctile(b_adj,2.5), prctile(b_adj,97.5), obs1{a}.cos_adj, prctile(b_cosadj,2.5), prctile(b_cosadj,97.5), ...
        obs1{a}.cos_startend, prctile(b_cosend,2.5), prctile(b_cosend,97.5));
    fprintf('  dissimilarity by lag 1..4: %s\n', mat2str(round(obs1{a}.lagrel,3)));
    fprintf('  DRIFT WITH LAG: Spearman rho = %.3f [%.3f, %.3f], slope = %.3f [%.3f, %.3f], epoch-order permutation p = %.4f (of 120 orders)\n', ...
        obs1{a}.rho, prctile(b_rho,2.5), prctile(b_rho,97.5), obs1{a}.slope, prctile(b_slope,2.5), prctile(b_slope,97.5), p_order);
    obs1{a}.p_order = p_order;
end

fprintf('\n--- Part 1: frontal minus auditory (neuron bootstrap, two-sided) ---\n');
names1 = {'reset','adj','rho','slope','cos_adj','cos_end'};
for i = 1:numel(names1)
    d = boot1(2).(names1{i}) - boot1(1).(names1{i});
    p = max(2*min(mean(d<=0), mean(d>=0)), 1/numel(d));
    fprintf('%-8s frontal - auditory = %+.3f  [%.3f, %.3f], p = %.4f\n', names1{i}, median(d), prctile(d,2.5), prctile(d,97.5), p);
end
d_slope = boot1(2).slope - boot1(1).slope;
p_slope = max(2*min(mean(d_slope<=0), mean(d_slope>=0)), 1/numel(d_slope));

%% ---- Step 3: Part 2 (letter C at positions 2/3/5; sequence-controlled) ----

fprintf('\n================ PART 2: letter C at positions 2/3/5, cross-sequence-group ================\n');
obs2 = struct(); boot2 = struct(); nullbar = cell(1,2);
for a = 1:2
    idx = area_idx{a};
    idx = idx(~isnan(DXn(idx,1)));
    N2(a) = numel(idx);
    DXbar = mean(mean(DXn(idx,:),2));
    DWbar = mean(mean(DWn(idx,:),2));
    Ebar  = mean(mean(SXn(idx,:),2));
    frac  = DXbar / DWbar;
    rel   = DXbar / Ebar;

    b_dx = nan(n_boot,1); b_dw = nan(n_boot,1); b_frac = nan(n_boot,1); b_rel = nan(n_boot,1);
    rowDX = mean(DXn(idx,:),2); rowDW = mean(DWn(idx,:),2); rowE = mean(SXn(idx,:),2);
    for b = 1:n_boot
        bi = randi(numel(idx), numel(idx), 1);
        b_dx(b) = mean(rowDX(bi)); b_dw(b) = mean(rowDW(bi));
        b_frac(b) = b_dx(b) / b_dw(b);
        b_rel(b)  = b_dx(b) / mean(rowE(bi));
    end
    nb = squeeze(mean(mean(DXnull(idx,:,:),3),1));    % 1 x n_perm, population null of mean DX
    nullbar{a} = double(nb(:));
    p_perm = (1 + sum(nullbar{a} >= DXbar)) / (1 + n_perm);

    obs2(a).DX = DXbar; obs2(a).DW = DWbar; obs2(a).frac = frac; obs2(a).p_perm = p_perm;
    boot2(a).dx = b_dx; boot2(a).dw = b_dw; boot2(a).frac = b_frac; boot2(a).rel = b_rel;

    fprintf('\n%s (%d neurons with >= %d trials in all 6 cells)\n', area_names{a}, N2(a), min_trials_c);
    for pr = 1:3
        fprintf('  position %d vs %d: cross-group distance = %.4f, within-group distance = %.4f\n', ...
            pos_c(pairs(pr,1)), pos_c(pairs(pr,2)), mean(DXn(idx,pr)), mean(DWn(idx,pr)));
    end
    fprintf('  CROSS-GROUP mean distance = %.4f [%.4f, %.4f] | permutation null median %.4f, 95th pct %.4f | permutation p = %.4f\n', ...
        DXbar, prctile(b_dx,2.5), prctile(b_dx,97.5), median(nullbar{a}), prctile(nullbar{a},95), p_perm);
    fprintf('  WITHIN-group mean distance = %.4f [%.4f, %.4f]\n', DWbar, prctile(b_dw,2.5), prctile(b_dw,97.5));
    fprintf('  fraction sequence-independent (cross / within) = %.3f [%.3f, %.3f]; cross-group distance relative to response energy = %.3f [%.3f, %.3f]\n', ...
        frac, prctile(b_frac,2.5), prctile(b_frac,97.5), rel, prctile(b_rel,2.5), prctile(b_rel,97.5));
end

fprintf('\n--- Part 2: frontal minus auditory ---\n');
d_dx = boot2(2).dx - boot2(1).dx;
p_dx = max(2*min(mean(d_dx<=0), mean(d_dx>=0)), 1/numel(d_dx));
nd   = nullbar{2} - nullbar{1};
p_dx_perm = (1 + sum(abs(nd) >= abs(obs2(2).DX - obs2(1).DX))) / (1 + n_perm);
fprintf('cross-group distance: frontal - auditory = %+.4f [%.4f, %.4f], bootstrap p = %.4f | permutation p = %.4f\n', ...
    median(d_dx), prctile(d_dx,2.5), prctile(d_dx,97.5), p_dx, p_dx_perm);
d_fr = boot2(2).frac - boot2(1).frac;
fprintf('fraction sequence-independent: frontal - auditory = %+.3f [%.3f, %.3f], p = %.4f\n', ...
    median(d_fr), prctile(d_fr,2.5), prctile(d_fr,97.5), max(2*min(mean(d_fr<=0), mean(d_fr>=0)), 1/numel(d_fr)));

%% ---- Step 4: pre-specified rule ----

crit_p1 = (obs1{2}.rho > 0) && (obs1{2}.p_order < 0.05) && (median(d_slope) > 0) && (p_slope < 0.05);
crit_p2 = (obs2(2).p_perm < 0.05) && (median(d_dx) > 0) && (p_dx < 0.05) && (p_dx_perm < 0.05);
fprintf('\n================ PRE-SPECIFIED READING ================\n');
fprintf('Part 1 (frontal progresses, auditory resets): frontal rho = %.3f (order-permutation p = %.4f); frontal-auditory slope p = %.4f  -> %s\n', ...
    obs1{2}.rho, obs1{2}.p_order, p_slope, ternary_str(crit_p1,'SUPPORTED','NOT SUPPORTED'));
fprintf('Part 2 (sequence-independent position signal, stronger in frontal): frontal permutation p = %.4f; frontal-auditory bootstrap p = %.4f, permutation p = %.4f  -> %s\n', ...
    obs2(2).p_perm, p_dx, p_dx_perm, ternary_str(crit_p2,'SUPPORTED','NOT SUPPORTED'));

%% ---- Step 5: plots (standard MATLAB plotting) ----

figure('Renderer','painters','Position',[100 100 1100 600]);
for a = 1:2
    subplot(2,3,a)
    imagesc(obs1{a}.cs, [0 1]); axis square; colorbar
    set(gca,'XTick',1:5,'YTick',1:5); xlabel('Epoch'); ylabel('Epoch');
    title([area_names{a} ': epoch cosine similarity'])
end
subplot(2,3,3); hold on
for a = 1:2
    lo = prctile(boot1(a).lag,2.5); hi = prctile(boot1(a).lag,97.5);
    patch([1:4 4:-1:1], [lo fliplr(hi)], cols_area(a,:), 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(1:4, obs1{a}.lagrel, '-o', 'Color', cols_area(a,:), 'LineWidth', 1.5);
end
set(gca,'XTick',1:4); xlabel('Epoch lag'); ylabel('Relative dissimilarity'); box off
title('Dissimilarity vs lag (95% CI, neuron bootstrap)')

subplot(2,3,4); hold on
for a = 1:2
    bar(a, obs2(a).DX, 0.4, 'FaceColor', cols_area(a,:));
    errorbar(a, obs2(a).DX, obs2(a).DX - prctile(boot2(a).dx,2.5), prctile(boot2(a).dx,97.5) - obs2(a).DX, 'k');
    plot(a + [-0.3 0.3], prctile(nullbar{a},95)*[1 1], 'k--');
end
set(gca,'XTick',1:2,'XTickLabel',area_names); ylabel('Cross-group distance'); box off
title('Part 2: sequence-independent (dashed = null 95th)')

subplot(2,3,5); hold on
for a = 1:2
    bar(a, obs2(a).DW, 0.4, 'FaceColor', cols_area(a,:));
    errorbar(a, obs2(a).DW, obs2(a).DW - prctile(boot2(a).dw,2.5), prctile(boot2(a).dw,97.5) - obs2(a).DW, 'k');
end
set(gca,'XTick',1:2,'XTickLabel',area_names); ylabel('Within-group distance'); box off
title('Part 2: including sequence context')

subplot(2,3,6); hold on
for a = 1:2
    bar(a, obs2(a).frac, 0.4, 'FaceColor', cols_area(a,:));
    errorbar(a, obs2(a).frac, obs2(a).frac - prctile(boot2(a).frac,2.5), prctile(boot2(a).frac,97.5) - obs2(a).frac, 'k');
end
set(gca,'XTick',1:2,'XTickLabel',area_names); ylabel('Cross / within'); box off
title('Fraction sequence-independent')


%% ================= Local functions =================

function m = traj_metrics(P)
% P: K x K population mean of cross-validated inner products between epochs.
K = size(P,1);
S = diag(P);
rel = nan(K,K); cs = nan(K,K);
for k = 1:K
    for j = 1:K
        if S(k) > 0 && S(j) > 0
            D = S(k) + S(j) - 2*P(k,j);
            rel(k,j) = D / (0.5*(S(k) + S(j)));
            cs(k,j)  = P(k,j) / sqrt(S(k)*S(j));
        end
    end
end
[pk, pj] = find(triu(true(K),1));
lag  = pj - pk;
relv = rel(sub2ind([K K], pk, pj));
m.rel = rel; m.cs = cs;
m.reset = mean(relv, 'omitnan');
m.adj   = mean(relv(lag==1), 'omitnan');
m.lagrel = nan(1,K-1);
for L = 1:K-1
    m.lagrel(L) = mean(relv(lag==L), 'omitnan');
end
ok = ~isnan(relv);
if sum(ok) >= 3
    m.rho = corr(lag(ok), relv(ok), 'Type','Spearman');
    pf = polyfit(lag(ok), relv(ok), 1);
    m.slope = pf(1);
else
    m.rho = nan; m.slope = nan;
end
m.cos_adj      = mean(cs(sub2ind([K K], pk(lag==1), pj(lag==1))), 'omitnan');
m.cos_startend = cs(1,K);
end


function s = ternary_str(cond, a, b)
if cond; s = a; else; s = b; end
end


function trial_bins = cell_feats(all_sdf, cond_value_all, cond_values, onset_ms, zero_col, time_win, ...
    bin_edges, n_bins, mu_fr, sd_fr)
trial_idx = ismember(cond_value_all, cond_values);
if isempty(cond_values) || ~any(trial_idx)
    trial_bins = [];
    return
end
w = all_sdf(trial_idx, zero_col + onset_ms + time_win);
raw = nan(size(w,1), n_bins);
for b = 1:n_bins
    raw(:,b) = mean(w(:, bin_edges(b)+1:bin_edges(b+1)), 2, 'omitnan');
end
raw = raw(all(~isnan(raw),2),:);
if isempty(raw)
    trial_bins = [];
else
    trial_bins = (raw - mu_fr) ./ sd_fr;
end
end
