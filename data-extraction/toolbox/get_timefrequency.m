function lfp_tf_out = get_timefrequency(lfp_aligned, ops)
%% Time frequency analysis
% Set up lfp_tf structure
lfp_tf.srate = 1000;
lfp_tf.pnts = length(ops.timewin);
lfp_tf.times = ops.timewin;

% wavelet parameters
min_freq = 2;
max_freq = 128;
num_frex = 30;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = -1:1/lfp_tf.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = lfp_tf.pnts;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 3; 


ch_extract = ops.ch_extract;

% Loop through each LFP channel
for lfp_ch_i = 1:length(ch_extract)
    lfp_ch = ch_extract(lfp_ch_i);

    % Clear loop variables
    lfp_tf.data = []; fft_data = [];

    % Call particular channel for LFP data
    lfp_tf.data = lfp_aligned.(['lfp_' int2str(lfp_ch)])(ops.tf_trials,:);

    % FFT of data (note: this doesn't change on frequency iteration)
    fft_data = fft(squeeze(lfp_tf.data(23,:,1)),n_conv_pow2);

    % Initialize output time-frequency data
    tf_data = zeros(length(frequencies),lfp_tf.pnts);

    % Run time-frequency power extraction
    for fi=1:length(frequencies)

        % create wavelet and get its FFT
        wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
        fft_wavelet = fft(wavelet,n_conv_pow2);

        % run convolution
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
        convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);

        % put power data into time-frequency matrix
        tf_data(fi,:) = abs(convolution_result_fft).^2;
    end

    % Make relevant baseline adjustments
    % define baseline period
    baselinetime = [ -1000 0 ]; % in ms

    % convert baseline window time to indices
    [~,baselineidx(1)]=min(abs(lfp_tf.times-baselinetime(1)));
    [~,baselineidx(2)]=min(abs(lfp_tf.times-baselinetime(2)));

    % dB-correct
    clear baseline_power dbconverted
    baseline_power = mean(tf_data(:,baselineidx(1):baselineidx(2)),2);
    dbconverted = 10*log10( bsxfun(@rdivide,tf_data,baseline_power));

    lfp_tf_power(:,:,lfp_ch_i) = dbconverted;
end


lfp_tf_out.power = lfp_tf_power;
lfp_tf_out.frequencies = frequencies;
lfp_tf_out.time = time;
