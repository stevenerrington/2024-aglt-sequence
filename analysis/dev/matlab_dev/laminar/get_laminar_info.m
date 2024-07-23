function laminar_info = get_laminar_info(electrode_ref, lfp_aligned, agl_t_event)

electrode_list = electrode_ref;
clear signal_out

%% Restructure data for spectrolaminar analysis

fprintf('Restructuring LFP data into the appropriate format   \n')
% - for each electrode
for electrode_i = 1:length(electrode_list)
    electrode_idx = electrode_list(electrode_i);
    n_trials = size(lfp_aligned.(['lfp_' int2str(electrode_i)]),1);

    % - across each trial
    for trial_i = 1:n_trials

        % - get the LFP data for the given
        signal_in = lfp_aligned.(['lfp_' int2str(electrode_idx)])(trial_i,:);

        % - save the relevant data in an output array for future use
        signal_out(electrode_i,:,trial_i) = signal_in; % nchans x trialtime x ntrials


    end
end

%% Find valid trials
valid_trial_idx = ~isnan(agl_t_event.stimulusOnset_ms);
signal_out_trl = signal_out(:,:,valid_trial_idx);

%% Get signal power (normalized and non-normalized; required FieldTrip)
fprintf('Calculating signal power   \n')

clear signal_normalized signal_nonnormalized
[signal_normalized, signal_nonnormalized] = lfp_to_normpower(signal_out_trl);


%% Run the main spectrolaminar extraction analysis code
fprintf('Running FLIP toolbox   \n')

% FLIP (Spectrolaminar) toolbox
[startinglowfreq,endinglowfreq,...
    startinghighfreq,endinghighfreq,...
    goodnessvalue,superficialchannel,deepchannel,...
    highfreqmaxchannel,lowfreqmaxchannel,...
    crossoverchannel] = ...
    FLIPAnalysis(signal_normalized,0:size(signal_normalized,1)-1,1:size(signal_normalized,2),1);

% Cross-contact correlation
fprintf('Running xcorr toolbox   \n')
xcontact_corr = D_LFPCORR_BASIC(signal_out_trl, 512, true);

% Current source density
fprintf('Running csd toolbox   \n')
csd_out = D_CSD_BASIC(signal_out_trl, 'spc', 0.2);

% Output data in a structure
laminar_info.lfp = signal_out_trl;


laminar_info.power.normalized = signal_normalized;
laminar_info.power.nonnormalized = signal_nonnormalized;

laminar_info.flip.goodnessvalue = goodnessvalue;
laminar_info.flip.superficialchannel = superficialchannel;
laminar_info.flip.deepchannel = deepchannel;
laminar_info.flip.highfreqmaxchannel = lowfreqmaxchannel;
laminar_info.flip.crossoverchannel = crossoverchannel;

laminar_info.xcontact = xcontact_corr;

laminar_info.csd = csd_out;


end
