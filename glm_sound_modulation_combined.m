%% Combined GLM Function: glm_sound_modulation_combined
%
% This version uses a combined model: firing_rate ~ sound + order_pos + exp_i
% and returns separate significance flags for sound, position, and interaction (optional).

function OUT = glm_sound_modulation_combined(sound_info_in, varargin)

p = inputParser;
p.addParameter('IncludeInteraction', false);
p.parse(varargin{:});
include_interaction = p.Results.IncludeInteraction;

%% --------------------- 1. Extract Data ----------------------------
reg_tbl = table;
sdf = cell2mat(sound_info_in(:,1));
reg_tbl.sound = string(sound_info_in(:,3));
reg_tbl.order_pos = string(sound_info_in(:,4));
reg_tbl.condition = string(sound_info_in(:,5));
reg_tbl.exp_i = (1:size(sdf,1))';
reg_tbl.valid = cell2mat(sound_info_in(:,9));

%% --------------------- 2. Windowed SDF ----------------------------
window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(sdf, window_size, window_shift);
[n_trials, n_times] = size(window_sdf);

analysis_win_idx = find(window_time >= 0 & window_time <= 413);

%% --------------------- 3. Init Arrays ------------------------------
% dynamic because # of coefficients depends on factor cardinality

% First get temporary model to read coefficient names
valid_trial_idx = find(reg_tbl.valid == 1);
tmp = reg_tbl(valid_trial_idx,:);
tmp.firing_rate = window_sdf(valid_trial_idx,1);

if include_interaction
    tmp_mdl = fitglm(tmp, 'firing_rate ~ sound*order_pos + exp_i');
else
    tmp_mdl = fitglm(tmp, 'firing_rate ~ sound + order_pos + exp_i');
end

coef_names = tmp_mdl.Coefficients.Properties.RowNames;
coef_count = numel(coef_names);

GLM_beta = nan(coef_count, n_times);
GLM_pval = nan(coef_count, n_times);
GLM_sig = false(coef_count, n_times);

%% --------------------- 4. Loop over Time Windows ------------------
for t = 1:n_times
    tbl = reg_tbl(valid_trial_idx,:);
    tbl.firing_rate = window_sdf(valid_trial_idx,t);

    if include_interaction
        mdl = fitglm(tbl, 'firing_rate ~ sound*order_pos + exp_i');
    else
        mdl = fitglm(tbl, 'firing_rate ~ sound + order_pos + exp_i');
    end

    GLM_beta(:,t) = mdl.Coefficients.Estimate;
    GLM_pval(:,t) = mdl.Coefficients.pValue;
end

%% --------------------- 5. FDR Across Time per Coefficient -----------
p_cutoff = 0.05;
for c = 1:coef_count
    q = mafdr(GLM_pval(c,:), 'BHFDR', true);
    GLM_sig(c,:) = q < p_cutoff;
end

%% --------------------- 6. Temporal Persistence Criterion -----------
signal_detect_length = 50; % ms
signal_detect_wins = signal_detect_length / window_shift;

coef_flag = false(coef_count,1);
for c = 1:coef_count
    [~, seg_len] = ZeroOnesCount(GLM_sig(c, analysis_win_idx));
    coef_flag(c) = any(seg_len >= signal_detect_wins);
end

%% --------------------- 7. Group Coefficients by Factor --------------
coef_sound = find(contains(coef_names, 'sound_'));
coef_pos   = find(contains(coef_names, 'order_pos_'));
coef_inter = find(contains(coef_names, ':'));

sound_flag = any(coef_flag(coef_sound));
pos_flag   = any(coef_flag(coef_pos));
if include_interaction
    inter_flag = any(coef_flag(coef_inter));
else
    inter_flag = [];
end

%% --------------------- 8. Package Outputs ---------------------------
OUT = struct();
OUT.GLM_beta = GLM_beta;
OUT.GLM_pval = GLM_pval;
OUT.GLM_sig  = GLM_sig;
OUT.coef_names = coef_names;
OUT.sound_flag = sound_flag;
OUT.pos_flag = pos_flag;
OUT.inter_flag = inter_flag;
OUT.window_sdf = window_sdf;
OUT.window_time = window_time;

end


%% =====================================================================
%% PERMUTATION TEST: Counts neurons responsive to sound/position
function OUT = permutation_test(all_units, varargin)

p = inputParser;
p.addParameter('Nperm', 1000);
p.addParameter('IncludeInteraction', false);
p.parse(varargin{:});
Nperm = p.Results.Nperm;
include_interaction = p.Results.IncludeInteraction;

Nunits = numel(all_units);

%% ---- Compute observed flags ----
obs_sound = false(Nunits,1);
obs_pos   = false(Nunits,1);
obs_inter = false(Nunits,1);

for u = 1:Nunits
    res = glm_sound_modulation_combined(all_units{u}, 'IncludeInteraction', include_interaction);
    obs_sound(u) = res.sound_flag;
    obs_pos(u)   = res.pos_flag;
    if include_interaction; obs_inter(u) = res.inter_flag; end
end

obs_count.sound = sum(obs_sound);
obs_count.pos   = sum(obs_pos);
if include_interaction; obs_count.inter = sum(obs_inter); end

%% ---- Null distributions ----
null_sound = zeros(Nperm,1);
null_pos   = zeros(Nperm,1);
null_inter = zeros(Nperm,1);

for p_i = 1:Nperm
    for u = 1:Nunits
        shuffled = all_units{u};
        idx = randperm(size(shuffled,1));
        shuffled(:,3) = shuffled(idx,3); % shuffle sound labels
        shuffled(:,4) = shuffled(idx,4); % shuffle positions

        res = glm_sound_modulation_combined(shuffled, 'IncludeInteraction', include_interaction);
        null_sound(p_i) = null_sound(p_i) + res.sound_flag;
        null_pos(p_i)   = null_pos(p_i) + res.pos_flag;
        if include_interaction
            null_inter(p_i) = null_inter(p_i) + res.inter_flag;
        end
    end
end

%% ---- p-values ----
pval.sound = mean(null_sound >= obs_count.sound);
pval.pos   = mean(null_pos >= obs_count.pos);
if include_interaction
    pval.inter = mean(null_inter >= obs_count.inter);
end

OUT = struct();
OUT.obs_count = obs_count;
OUT.null_sound = null_sound;
OUT.null_pos   = null_pos;
OUT.null_inter = null_inter;
OUT.pval = pval;

end


%% =====================================================================
%% VISUALIZATION CODE
function plot_population_results(obs_count, null_dist)
figure; 
subplot(1,2,1);
histogram(null_dist);
hold on;
xline(obs_count, 'LineWidth',2);
title('Null distribution');
xlabel('# units'); ylabel('count');
subplot(1,2,2);
bar([obs_count]); ylabel('# units'); title('Observed');
end