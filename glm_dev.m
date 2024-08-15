
clear glm_output encoding_flag encoding_beta

parfor neuron_i = 1:size(spike_log,1)
    fprintf('Neuron %i of %i \n', neuron_i, size(spike_log,1))

    [glm_output{neuron_i}, encoding_flag(neuron_i,:), encoding_beta(neuron_i,:)] =...
        glm_sound_modulation(neuron_sdfsound_out{neuron_i});
end



neurons_in = find(sum(encoding_flag(:,[1:4]),2) > 0);


for neuron_i = 1:10
figuren;
hold on
plot(ops.sound_sdf_window,norm_fr_soundA(neurons_in(neuron_i),:))
plot(ops.sound_sdf_window,norm_fr_soundC(neurons_in(neuron_i),:))
plot(ops.sound_sdf_window,norm_fr_soundD(neurons_in(neuron_i),:))
plot(ops.sound_sdf_window,norm_fr_soundF(neurons_in(neuron_i),:))
plot(ops.sound_sdf_window,norm_fr_soundG(neurons_in(neuron_i),:))
title(neuron_i)
    
end

figuren;
subplot(1,2,1); hold on
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundA(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundC(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundD(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundF(intersect(neurons_in,frontal_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundG(intersect(neurons_in,frontal_neuron_idx),:)),50))

subplot(1,2,2); hold on
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundA(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundC(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundD(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundF(intersect(neurons_in,auditory_neuron_idx),:)),50))
plot(ops.sound_sdf_window,smooth(nanmean(norm_fr_soundG(intersect(neurons_in,auditory_neuron_idx),:)),50))

length(intersect(neurons_in,frontal_neuron_idx))
length(intersect(neurons_in,auditory_neuron_idx))

%%
neuron_i = 1499;

clear glm_output encoding_flag encoding_beta

 [glm_output, encoding_flag, encoding_beta] =...
        glm_sound_modulation(neuron_sdfsound_out{neuron_i})