nNeurons = size(spike_log,1);

%%
clear results
parfor neuron_i = 1:nNeurons
    sound_info_in = sdf_soundAlign_data{neuron_i}

    %% Extract: get relevant data for GLM table
    reg_tbl = table;

    sound_info_in_sdf = cell2mat(sound_info_in(:,1));

    reg_tbl.identity = string(sound_info_in(:,3));
    reg_tbl.position = string(sound_info_in(:,4));
    reg_tbl.condition = string(sound_info_in(:,5));
    reg_tbl.exp_i = [1:size(sound_info_in_sdf,1)]';
    reg_tbl.valid = cell2mat(sound_info_in(:,9));


    %% Setup spike data into GLM
    baseline_idx = 200+(-100:0);
    mean_fr = nanmean(nanmean(sound_info_in_sdf(:,baseline_idx),2));
    std_fr = nanstd(nanmean(sound_info_in_sdf(:,baseline_idx),2));
    z_sdf = (sound_info_in_sdf - mean_fr) ./ std_fr;

    window_size = 100;
    window_shift = 10;
    [window_sdf, window_time] = movaverage_sdf(sound_info_in_sdf, window_size, window_shift);
    % 2023-08-26, 23h25: I tested this and it does the job correctly. The
    % z-scored profile looks exactly like the raw profile, and the movaverage
    % sdf looks like a smoothed version of the raw/zscored sdf.

    [n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

    baseline_win_idx = find(window_time >= -200 & window_time <= 0); % Find the relevant indicies for the timepoints of interest
    analysis_win_idx = find(window_time >= 0 & window_time <= 413); % Find the relevant indicies for the timepoints of interest
    win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

    reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

    %%
    valid_trials = [];
    valid_trials = find(reg_tbl.valid == 1 & ~strcmp(reg_tbl.position,'position_0') & strcmp(reg_tbl.condition,'nonviol'));

    reg_tbl = reg_tbl(valid_trials,:);
    reg_tbl.identity = categorical(reg_tbl.identity);
    reg_tbl.position = categorical(reg_tbl.position);

    %%
    mdl = fitglm(reg_tbl, 'win_fr ~ identity + position + identity:position','Distribution','Poisson');
    % Extract p-values
    p_identity   = mdl.Coefficients.pValue(contains(mdl.CoefficientNames,'identity') & ~contains(mdl.CoefficientNames,'position'));
    p_position   = mdl.Coefficients.pValue(contains(mdl.CoefficientNames,'position') & ~contains(mdl.CoefficientNames,'identity'));
    p_interaction= mdl.Coefficients.pValue(contains(mdl.CoefficientNames,'identity') & contains(mdl.CoefficientNames,'position'));

    % Combine p-values for overall effects in each term
    p_identity      = min(p_identity);
    p_position      = min(p_position);
    p_interaction   = min(p_interaction);

    % Save
    results(neuron_i).p_identity    = p_identity;
    results(neuron_i).p_position    = p_position;
    results(neuron_i).p_interaction = p_interaction;

end

%%
alpha = 0.05;

for n = 1:nNeurons
    pi = results(n).p_identity;
    pp = results(n).p_position;
    px = results(n).p_interaction;

    if (pi < alpha) && (pp >= alpha)
        results(n).type = 'Identity';
        
    elseif (pp < alpha) && (pi >= alpha)
        results(n).type = 'Position';
        
    elseif (pi < alpha) && (pp < alpha) && (px >= alpha)
        results(n).type = 'Both (additive)';
        
    elseif (px < alpha)
        results(n).type = 'Mixed-selectivity (interaction)';
        
    else
        results(n).type = 'Unselective';
    end
end

%%
results = struct2table(results);

sum(strcmp(results.type(glm_sig_units),'Unselective'))
sum(strcmp(results.type(glm_sig_units),'Identity'))
sum(strcmp(results.type(glm_sig_units),'Position'))
sum(strcmp(results.type(glm_sig_units),'Both (additive)'))
sum(strcmp(results.type(glm_sig_units),'Mixed-selectivity (interaction)'))

