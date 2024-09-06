% Define a list of predefined sequences
seq_list{1} = {'Start','A', 'C', 'F', 'C','X'};
seq_list{2} = {'Start','A', 'C', 'F', 'C','G','X'};
seq_list{3} = {'Start','A', 'C', 'G', 'F'};
seq_list{4} = {'Start','A', 'C', 'G', 'F','C','G','X'};
seq_list{5} = {'Start','A', 'D', 'C', 'F','X'};
seq_list{6} = {'Start','A', 'D', 'C', 'F','C','X'};
seq_list{7} = {'Start','A', 'D', 'C', 'F','C','G','X'};
seq_list{8} = {'Start','A', 'D', 'C', 'G','F','C','X'};

% Set the random seed for reproducibility of the randomization
rng(1, 'twister')

% Define the total number of repetitions needed
reps = 160;

% Create an index list for the 8 sequences
seq_indx = 1:8;

% Initialize an empty array to store the randomized sequence indices
trial_order = [];

% Loop until 160 elements (trials) are collected in 'trial_order'
while length(trial_order) < reps
    % Randomly shuffle the sequence indices and append to 'trial_order'
    trial_order = [trial_order, seq_indx(randperm(length(seq_indx)))];
end

% Define the elements used in sequences
elements = {'Start','A','C','D','F','G','X'};

% Initialize a transitional probability matrix (7x7 for 7 elements)
transitional_prob_array = zeros(7,7);

% Initialize a positional probability matrix (7x7 for 7 elements)
positional_prob_array = zeros(7,8);

% Loop over all trials to count transitions between elements
for trial_i = 1:reps
    % Get the sequence corresponding to the current trial
    seq_in = seq_list{trial_order(trial_i)};
    
    % Loop through elements in the sequence and track transitions
    for ele_i = 2:length(seq_in)
        % Find the index of the previous element in 'elements'
        last_ele_index = find(strcmp(elements, seq_in{ele_i-1}));
        
        % Find the index of the current element in 'elements'
        next_ele_index = find(strcmp(elements, seq_in{ele_i}));
        
        % Increment the transition count between the previous and current element
        transitional_prob_array(last_ele_index, next_ele_index) = ...
            transitional_prob_array(last_ele_index, next_ele_index) + 1;

        % Increment the positional count
        positional_prob_array(next_ele_index, ele_i) = ...
            positional_prob_array(next_ele_index, ele_i) + 1;
    end
end

% Create color scales for backward and forward heatmaps
colorscale_backward = cbrewer('seq', 'Greens', 10);  % Green color scale for backward transitions
colorscale_forward = cbrewer('seq', 'Reds', 10);    % Red color scale for forward transitions
colorscale_position = cbrewer('seq', 'Purples', 10);    % Grey color scale for positional probabilities

% Create a figure with white background
figure('Renderer', 'painters', 'Position', [100 100 2000 400], 'Color', 'W');

% Plot the backward transitional probability heatmap
subplot(1, 3, 1)
backward_heatmap = heatmap(round((transitional_prob_array ./ sum(transitional_prob_array, 1))*100,1));
backward_heatmap.XDisplayLabels = elements;  % Set custom X-axis labels
backward_heatmap.YDisplayLabels = elements;  % Set custom Y-axis labels
xlabel('Element e')                          % X-axis label
ylabel('Element e-1')                        % Y-axis label
title('Backward transitional probability')    % Title of the plot
grid off                                     % Disable grid lines
colormap(backward_heatmap, colorscale_backward);  % Apply green colormap

% Plot the forward transitional probability heatmap
subplot(1, 3, 2)
forward_heatmap = heatmap(round((transitional_prob_array ./ sum(transitional_prob_array, 2))*100,1));
forward_heatmap.XDisplayLabels = elements;   % Set custom X-axis labels
forward_heatmap.YDisplayLabels = elements;   % Set custom Y-axis labels
xlabel('Element e+1')                        % X-axis label
ylabel('Element e')                          % Y-axis label
title('Forward transitional probability')    % Title of the plot
grid off                                     % Disable grid lines
colormap(forward_heatmap, colorscale_forward);  % Apply red colormap

% Plot the positional probability heatmap
subplot(1, 3, 3)
positional_heatmap = heatmap(round((positional_prob_array ./ sum(positional_prob_array, 2))*100,1));
positional_heatmap.XDisplayLabels = {'Start', '1st','2nd','3rd','4th','5th','6th','7th'};   % Set custom X-axis labels
positional_heatmap.YDisplayLabels = elements;   % Set custom Y-axis labels
xlabel('Position')                           % X-axis label
ylabel('Element')                            % Y-axis label
title('Positional probability')              % Title of the plot
grid off                                     % Disable grid lines
colormap(positional_heatmap, colorscale_position);  % Apply grey colormap