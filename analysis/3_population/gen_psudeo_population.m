function subsamples = gen_psudeo_population(pop_index, n_obs_psuedopop)

% Generate psudeo populations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_obs_pop = size(pop_index,1);

rng(1,"twister");
% Initialize the list of numbers
n_obs_list = pop_index;


% Define the number of subsamples you want
x = floor(n_obs_pop/n_obs_psuedopop);

% Preallocate a cell array to store the subsamples
subsamples = cell(1, x);

% Loop to create x subsamples
for i = 1:x
    % Check if there are enough remaining numbers to sample
    if length(n_obs_list) < n_obs_psuedopop
        error('Not enough numbers left to create the next subsample.');
    end
    
    % Sample without replacement
    sampled_values = randsample(n_obs_list, n_obs_psuedopop);
    
    % Store the sampled values
    subsamples{i} = sampled_values;
    
    % Remove the sampled values from the list
    n_obs_list = setdiff(n_obs_list, sampled_values);
end