










%% Analysis: curate neurons with significant modulation
for neuron_i = 1:length(sig_flag)
    if sum(sig_flag(neuron_i,:)) == 0
        type_flag{neuron_i,1} = 'nonsig';
    else
        if all(dir_flag(sig_flag(neuron_i,:)) == 0)
            type_flag{neuron_i,1} = 'neg';
        elseif all(dir_flag(sig_flag(neuron_i,:)) == 1)
            type_flag{neuron_i,1} = 'pos';
        else
            type_flag{neuron_i,1} = 'mix';
        end
    end
end

sig_units = find(sum(sig_flag,2) > 0);

sig_pos_units = find(strcmp(type_flag,'pos'));
sig_neg_units = find(strcmp(type_flag,'neg'));
sig_mix_units = find(strcmp(type_flag,'mix'));
nonsig_units = find(strcmp(type_flag,'nonsig'));

sig_pos_units_auditory = intersect(sig_pos_units,auditory_neuron_idx);
sig_pos_units_frontal = intersect(sig_pos_units,frontal_neuron_idx);
sig_neg_units_auditory = intersect(sig_neg_units,auditory_neuron_idx);
sig_neg_units_frontal = intersect(sig_neg_units,frontal_neuron_idx);

sig_units_auditory = intersect(sig_units,auditory_neuron_idx);
sig_units_frontal = intersect(sig_units,frontal_neuron_idx);






