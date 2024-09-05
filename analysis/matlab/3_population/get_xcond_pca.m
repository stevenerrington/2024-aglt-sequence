function [pc_out, pc_shuf_out] = get_xcond_pca(signal_in)

n_conditions = length(signal_in);
rng(1,"twister");

% Concatenate signals across conditions %%%%%%%%%%%%%%%%%%%%%%%%%%
signal_in_conc = [];
signal_in_shuffled_conc = [];

for cond_i = 1:n_conditions
    signal_in_conc = [signal_in_conc, signal_in{cond_i}];

    shuffle_temp = [];
    for neuron_i = 1:size(signal_in{cond_i},1)
        shuffle_temp(neuron_i,:) = signal_in{cond_i}(neuron_i,randperm(size(signal_in{cond_i},2)));
    end

    signal_in_shuffled_conc = [signal_in_shuffled_conc, shuffle_temp];
end

% Parse condition indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond_label = [];
for cond_i = 1:n_conditions
    cond_label = [cond_label, repmat(cond_i,1,size(signal_in{cond_i},2))];
end

for cond_i = 1:n_conditions
    cond_idx{cond_i} = find(cond_label == cond_i);
end

% Perform PCA on regular and shuffled data %%%%%%%%%%%%%%%%%%%%%%%%
[coeff, pcs, latent, ~, var_exp, ~] = pca(signal_in_conc');
[coeff_shuffled, pcs_shuffled, latent_shuffled, ~, var_exp_shuffled, ~] = pca(signal_in_shuffled_conc');

for cond_i = 1:n_conditions
    pc_out{cond_i} = pcs(cond_idx{cond_i},:)';
    pc_shuf_out{cond_i} = pcs_shuffled(cond_idx{cond_i},:)';
end