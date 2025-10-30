function [f, coh_matrix] = compute_probe_coherence(Data, fs, win)
    % Data: [32 x nTimes x nTrials]
    % fs: sampling rate
    % win: window length (samples) for mscohere

    nCh1 = 16;
    nCh2 = 16;
    nTrials = size(Data, 3);

    % preallocate
    [Cxy, f] = mscohere(squeeze(Data(1,:,1)), squeeze(Data(17,:,1)), ...
                        hanning(win), win/2, win, fs);
    coh_matrix = zeros(nCh1, nCh2, length(Cxy));

    for i = 1:nCh1
        for j = 1:nCh2
            trial_coh = zeros(nTrials, length(Cxy));
            for tr = 1:nTrials
                sig1 = squeeze(Data(i, :, tr));
                sig2 = squeeze(Data(j+16, :, tr));
                [Cxy, f] = mscohere(sig1, sig2, hanning(win), win/2, win, fs);
                trial_coh(tr, :) = Cxy;
            end
            coh_matrix(i, j, :) = mean(trial_coh, 1);
        end
    end
end
