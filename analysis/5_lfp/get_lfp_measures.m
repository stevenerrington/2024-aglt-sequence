function [phase_analysis_out] = get_lfp_measures(lfp_data_in, alignTimes, filters, timewin)
% This function processes local field potential (LFP) data to extract phase and amplitude measures.
% It applies bandpass filtering, extracts instantaneous phase using Hilbert transform,
% aligns data to specified event times, and computes phase-locking values (PLV).

% Inputs:
%   lfp_data_in  - Matrix of LFP signals (channels x time)
%   alignTimes   - Time points to align trials around (e.g., stimulus onset times)
%   filters      - Cell array of bandpass filters for different frequency bands
%   timewin      - Time window for trial extraction around alignTimes

% Output:
%   phase_analysis_out - Struct containing phase, amplitude, and phase-locking values

n_channels = size(lfp_data_in, 1); % Number of recording channels
filter_labels = {'theta','alpha','beta','gamma'};

% Loop over each bandpass filter (corresponding to different frequency bands)
for fq = 1:numel(filters)
    for ch_i = 1:n_channels
        
        % Apply the bandpass filter to the LFP signal for the current channel
        fq_filt = filters{fq}; % Get filter coefficients for the current frequency band
        fq_signal = fftfilt(fq_filt.Coefficients, double(lfp_data_in(ch_i, :))); % Filter the LFP signal
        
        % Align the filtered LFP signal to the given event times
        [freq_lfp] = align_lfp_trials(fq_signal, alignTimes, timewin);
        
        % Adjust time to account for filter lag (this can be important for accurate phase estimation)
        AD_delay = 0; % Set to zero if not correcting for acquisition delays
        LFPtime = 1:numel(lfp_data_in(ch_i, :)); % Define time vector (in ms)
        lag = 1 + (1 + 1) * AD_delay / 2; % Compute lag adjustment (based on acquisition delay if present)
        adjusted_sig = fq_signal(lag:end); % Adjusted signal after correcting for lag
        
        % Extract the instantaneous phase using Hilbert Transform
        fq_phase = angle(hilbert(adjusted_sig)); % Compute phase in radians (-pi to pi)
        
        % Align phase data to the event times
        [phase_trials] = align_lfp_trials(fq_phase, alignTimes, timewin);
        
        % Compute phase alignment across trials (e.g., Phase Locking Value, PLV)
        % Remove trials with NaNs before computing alignment
        [phase_alignment] = compute_LFP_phase_alignment(phase_trials(~isnan(phase_trials(:, 1)), :));
        
        % Extract the amplitude envelope of the signal (magnitude of analytic signal)
        fq_amp = abs(hilbert(adjusted_sig));
        smoothed_amp = smooth(fq_amp, 50); % Smooth amplitude to reduce noise
        
        % Align amplitude data to event times
        [lfp_amp] = align_lfp_trials(smoothed_amp, alignTimes, timewin);
        
        % Store results in structured output for this frequency band and channel
        phase_analysis_out.([filter_labels{fq} '_lfp']){ch_i, 1} = freq_lfp; % Filtered LFP
        phase_analysis_out.([filter_labels{fq} '_phase']){ch_i, 1} = phase_trials; % Phase values
        phase_analysis_out.([filter_labels{fq} '_plv']){ch_i, 1} = phase_alignment; % Phase alignment (PLV)
        phase_analysis_out.([filter_labels{fq} '_amplitude']){ch_i, 1} = lfp_amp; % Amplitude envelope
        
    end % End channel loop
end % End frequency loop
