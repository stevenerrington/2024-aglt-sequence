ppc_out_session = {};
plv_out_session = {};

for session_i = 1:size(ephysLog,1)
    session = ephysLog.session{session_i};
    fprintf('Session %i of %i \n', session_i, size(ephysLog,1));
    % Load the corresponding session data from the specified directory
    session_data = load(fullfile(dirs.mat_data, session));

    session_entries = []; session_entries = find(strcmp(session,spike_log.session));
    session_spike_log = spike_log(session_entries,:);
    session_aud_neurons = []; session_aud_neurons = find(strcmp(session_spike_log.area,'auditory'));
    session_frontal_neurons = []; session_frontal_neurons = find(strcmp(session_spike_log.area,'frontal'));

    Fs = 1000;

    probe1_area = ephysLog.area_label_sec{session_i};

    switch probe1_area
        case 'auditory'
            lfp_range = [1:16];
        case 'frontal'
            lfp_range = [17:32];
    end

    for lfp_ch_idx = 1:16

        lfp_ch = lfp_range(lfp_ch_idx);

        lfp = []; lfp = session_data.lfp(lfp_ch,:);

        f1 = 4; f2 = 9;
        [b,a] = butter(3,[f1 f2]/(Fs/2));   % 3rd-order Butterworth
        lfp_filt = []; lfp_filt = filtfilt(b,a,lfp);       % zero-phase filtering

        lfp_phase = []; lfp_phase = angle(hilbert(lfp_filt));  % radians, same length as LFP


        for neuron_i = 1:length(session_frontal_neurons)
            spikes = []; spikes = round(session_data.spikes.time.(session_spike_log.unitDSP{session_frontal_neurons(neuron_i)}))/1000;

            % Convert spike times to sample indices
            spike_idx = round(spikes*Fs);
            spike_idx(spike_idx<1 | spike_idx>length(lfp_phase)) = [];

            % Extract phase at spike times
            spike_phases = []; spike_phases = lfp_phase(spike_idx);

            mvl(lfp_ch_idx, neuron_i) = abs(mean(exp(1i*spike_phases)));

            % Pairwise Phase Consistency
            nSpikes = length(spike_phases);
            ppc = 0;
            for k = 1:nSpikes-1
                for l = k+1:nSpikes
                    ppc = ppc + cos(spike_phases(k)-spike_phases(l));
                end
            end
            ppc = (2/(nSpikes*(nSpikes-1))) * ppc;

            ppc_out(lfp_ch_idx, neuron_i) = ppc;

        end

    end

    plv_out_session{session_i} = mvl;
    ppc_out_session{session_i} = ppc_out;

end


figure;
plv_out_session{session_i}