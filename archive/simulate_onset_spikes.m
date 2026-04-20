function [spike_times, unit_ids, params] = simulate_onset_spikes(n_units, n_trials)
% SIMULATE_ONSET_SPIKES Generate biologically plausible spike trains
% driven by auditory element onsets
%
% Inputs:
%   n_units  - Number of simulated neurons (default: 10)
%   n_trials - Number of repetitions of the sequence (default: 20)
%
% Outputs:
%   spike_times - Cell array of spike times for each unit
%   unit_ids    - Unit assignments for pooled spike analysis
%   params      - Structure containing stimulus and simulation parameters

if nargin < 1, n_units = 10; end
if nargin < 2, n_trials = 20; end

%% Stimulus parameters
params.element_duration = 0.413;  % seconds
params.isi = 0.150;               % inter-stimulus interval (seconds)
params.n_elements = 5;            % number of elements in sequence
params.n_trials = n_trials;

% Calculate element onset times
element_period = params.element_duration + params.isi;
params.onset_times = (0:params.n_elements-1) * element_period;

%% Neural response parameters
params.response_latency_mean = 0.015;  % 15ms mean latency
params.response_latency_std = 0.005;   % 5ms std of latency
params.response_window = 0.100;        % 100ms response window after onset
params.refractory_period = 0.002;      % 2ms absolute refractory period

%% Generate spikes for each unit
spike_times = cell(n_units, 1);
rng(42); % For reproducibility

for unit = 1:n_units
    % Unit-specific parameters (heterogeneity across population)
    unit_params = generate_unit_params(unit, n_units);
    
    % Generate spikes for this unit across all trials
    unit_spikes = [];
    
    for trial = 1:n_trials
        trial_offset = (trial - 1) * (params.onset_times(end) + 1.0);
        
        for elem = 1:params.n_elements
            % Adaptation: reduce response for later elements
            adaptation_factor = unit_params.adaptation_rate ^ (elem - 1);
            
            % Current element onset time
            onset_time = trial_offset + params.onset_times(elem);
            
            % Generate response latency for this trial/element
            latency = params.response_latency_mean + ...
                      params.response_latency_std * randn();
            latency = max(latency, 0.001); % Ensure positive
            
            % Response timing parameters
            response_start = onset_time + latency;
            response_end = response_start + params.response_window;
            
            % Generate spikes in response window
            element_spikes = generate_response_spikes(...
                response_start, response_end, ...
                unit_params.peak_rate * adaptation_factor, ...
                unit_params.temporal_profile, ...
                params.refractory_period);
            
            unit_spikes = [unit_spikes; element_spikes];
        end
    end
    
    % Sort spike times
    spike_times{unit} = sort(unit_spikes);
end

%% Create pooled spike data with unit IDs (useful for raster plots)
unit_ids = [];
all_spikes = [];
for unit = 1:n_units
    unit_ids = [unit_ids; unit * ones(length(spike_times{unit}), 1)];
    all_spikes = [all_spikes; spike_times{unit}];
end

%% Visualization
visualize_responses(spike_times, params);

end

%% Helper function: Generate unit-specific parameters
function unit_params = generate_unit_params(unit_id, n_units)
    % Create diversity in neural responses
    rng(unit_id); % Different for each unit, but reproducible
    
    % Peak firing rate (spikes/sec) - varies across units
    unit_params.peak_rate = 20 + 60 * rand();
    
    % Adaptation rate (0-1, where 1 = no adaptation)
    unit_params.adaptation_rate = 0.7 + 0.3 * rand();
    
    % Temporal profile: 'phasic' (transient) or 'tonic' (sustained)
    unit_params.temporal_profile = rand() > 0.5;
end

%% Helper function: Generate spikes for a single response
function spikes = generate_response_spikes(t_start, t_end, peak_rate, ...
                                           is_phasic, refrac_period)
    spikes = [];
    
    % Time step for simulation
    dt = 0.001; % 1ms
    t = t_start:dt:t_end;
    
    % Create firing rate profile over time
    if is_phasic
        % Phasic: strong onset response that decays
        tau_decay = 0.020; % 20ms decay time constant
        rate_profile = peak_rate * exp(-(t - t_start) / tau_decay);
    else
        % Tonic: sustained response with slight decay
        tau_decay = 0.100; % 100ms decay time constant
        rate_profile = peak_rate * exp(-(t - t_start) / tau_decay);
    end
    
    % Generate spikes using inhomogeneous Poisson process
    last_spike_time = -inf;
    
    for i = 1:length(t)
        % Check refractory period
        if (t(i) - last_spike_time) < refrac_period
            continue;
        end
        
        % Spike probability in this time bin
        spike_prob = rate_profile(i) * dt;
        
        if rand() < spike_prob
            spikes = [spikes; t(i)];
            last_spike_time = t(i);
        end
    end
end

%% Visualization function
function visualize_responses(spike_times, params)
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot 1: Raster plot for first 3 trials
    subplot(3, 1, 1);
    n_units = length(spike_times);
    trial_duration = params.onset_times(end) + 1.0;
    
    for unit = 1:n_units
        spikes = spike_times{unit};
        % Only plot first 3 trials
        trial_spikes = spikes(spikes < 3 * trial_duration);
        
        for spike = trial_spikes'
            line([spike, spike], [unit-0.4, unit+0.4], 'Color', 'k', 'LineWidth', 1);
        end
    end
    
    % Mark element onsets
    for trial = 0:2
        for onset = params.onset_times
            xline(trial * trial_duration + onset, 'r--', 'Alpha', 0.3);
        end
    end
    
    ylim([0, n_units+1]);
    xlim([0, 3 * trial_duration]);
    ylabel('Unit #');
    title('Raster Plot (First 3 Trials)');
    grid on;
    
    % Plot 2: PSTH aligned to first element onset
    subplot(3, 1, 2);
    bin_size = 0.005; % 5ms bins
    bins = -0.05:bin_size:0.3;
    psth = zeros(size(bins));
    
    for unit = 1:n_units
        spikes = spike_times{unit};
        for trial = 0:(params.n_trials-1)
            trial_offset = trial * trial_duration;
            aligned_spikes = spikes - (trial_offset + params.onset_times(1));
            psth = psth + histcounts(aligned_spikes, [bins, bins(end)+bin_size]);
        end
    end
    
    % Convert to firing rate
    psth = psth / (n_units * params.n_trials * bin_size);
    
    bar(bins, psth, 'k');
    xline(0, 'r--', 'LineWidth', 2);
    xlabel('Time from stimulus onset (s)');
    ylabel('Firing rate (spikes/s)');
    title('PSTH aligned to first element');
    grid on;
    
    % Plot 3: Response to each element (adaptation)
    subplot(3, 1, 3);
    element_responses = zeros(1, params.n_elements);
    
    for elem = 1:params.n_elements
        response_window = [0.010, 0.060]; % 10-60ms post-onset
        
        for unit = 1:n_units
            spikes = spike_times{unit};
            for trial = 0:(params.n_trials-1)
                trial_offset = trial * trial_duration;
                onset = trial_offset + params.onset_times(elem);
                
                count = sum(spikes >= onset + response_window(1) & ...
                           spikes < onset + response_window(2));
                element_responses(elem) = element_responses(elem) + count;
            end
        end
    end
    
    % Normalize
    element_responses = element_responses / (n_units * params.n_trials);
    
    bar(1:params.n_elements, element_responses, 'FaceColor', [0.3, 0.6, 0.9]);
    xlabel('Element number');
    ylabel('Mean spike count (10-60ms window)');
    title('Adaptation across sequence elements');
    grid on;
    
    sgtitle('Simulated Neural Responses to Auditory Sequence');
end
