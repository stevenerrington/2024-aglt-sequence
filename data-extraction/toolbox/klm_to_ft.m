function ft_data = klm_to_ft(lfp_data, ops)

% Create the FieldTrip data structure
ft_data = [];
ft_data.label = arrayfun(@(x) ['lfp_' num2str(x)], 1:size(lfp_data, 1), 'UniformOutput', false)'; % Create channel labels
ft_data.fsample = 1000; % Sampling frequency in Hz
ft_data.time = arrayfun(@(x) ops.timewin/1000, 1:size(lfp_data, 3), 'UniformOutput', false); % Time axis
ft_data.trial = arrayfun(@(x) squeeze(lfp_data(:, :, x)), 1:size(lfp_data, 3), 'UniformOutput', false); % Trials

% Add additional information if needed (e.g., trialinfo)
ft_data.trialinfo = ones(size(lfp_data, 3), 1); % Example trial info, all ones

end