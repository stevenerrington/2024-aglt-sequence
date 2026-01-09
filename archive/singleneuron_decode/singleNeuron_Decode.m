
%% Run identity & position decoding for all neurons

nNeurons = numel(sdf_soundAlign_data);

acc_identity = nan(nNeurons,1);
acc_position = nan(nNeurons,1);

clear decode_details

parfor n = 1:nNeurons
    fprintf('Neuron %i of %i\n', n, nNeurons);

    [acc_identity(n), acc_position(n), decode_details(n)] = ...
        decode_identity_position_singleNeuron(sdf_soundAlign_data{n});
end


%%
for n = 1:nNeurons
    singleNeuron_decode_identity(n,:) = decode_details(n).identity.classAcc';
    singleNeuron_decode_position(n,:) = decode_details(n).position.classAcc';
end





singleNeuron_decode_identity(neuron_class.frontal.all,:)
singleNeuron_decode_position(neuron_class.frontal.all,:)