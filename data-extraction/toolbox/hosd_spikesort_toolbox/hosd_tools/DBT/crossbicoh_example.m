% A demonstration of univariate, bivariate and trivariate bicoherence using
% pspect2
%
% C. Kovach 2021

N = 1e5; %Test signal duration in samples
Fs = 1e3; %Test signal sampling rate.
f_slow = [4 7]; % Band for the slow oscillation (SO)
f_fast = [80 120]; %Band for the fast oscillation (FO)
bandwidth = 2; % Bandwidth for the analysis window of the slow oscillation
bandwidth_ratio=1; %Bandwidth for the fast oscillation is bandwidth*bandwidth_ratio. 
                   %Note that setting this greater than 1 is not necessary,
                   %but doing so smooths along the FO dimension in a way that
                   %might improve sensitivity to PAC while reducing the size
                   %of the bicoherence array along the FO
                   %dimension.
slow_lowpass = 150; %Analyze slow oscillations only up to this frequency. 
                    %This is kept large for the purpose of the demonstration,
                    %To illustrate symmetry properties.
                    
noisamp = 1; %Noise amplitude

% simulate = 'univariate';
% simulate = 'bivariate';
% simulate = 'trivariate';

%%% Cases to simulate
simulate = {'univariate','bivariate','trivariate'};


if ~iscell(simulate)
    simulate={simulate};
end

%%% Slow oscillation
s1 = hilbert(zscore(bpfilt(randn(N,1),[Fs f_slow])));
% s1b = hilbert(zscore(bpfilt(randn(N,1),[Fs f_slow+5])));
%%% Fast oscillation
s2 = hilbert(zscore(bpfilt(randn(N,1),[Fs f_fast])));

%%% s3.*conj(s2) is modulated by s1
s3 = conj(s1).*s2*exp(1i*2*pi*rand)./(abs(s2).^2);
% s3 = s3+s1b.*s2*exp(1i*2*pi*rand)./(abs(s2).^2);
% s1=s1+s1b;

%%% Now create two better behaved signals, s4 and s5, from s2 and s3 with a 
%%% common envelope that is the geometric mean of the respective envelopes. 
%%% The phase difference between s4 and s5 is modulated by s1 while the power
%%% envelope of s4+s5 is likewise modulated.
s4 = zscore(s3.*sqrt(abs(s2)./abs(s3)));
s5 = zscore(s2.*sqrt(abs(s3)./abs(s2)));

figure('units','normalized','position',[0.05    0.5   0.9    0.4])

for k = 1:length(simulate)
    switch simulate{k}
        case {1,'univariate'}

            %%%Univariate signal composed of s1+s4+s5 and noise.
            %%%Univariate bicoherence for a real-valued signal has a 12-fold
            %%%symmetry under  the 6 permutations of all 3 frequency bands 
            %%% (f1,f2,-f1-f2) and sign reversal of all frequencies.
            x1 = real(s1+ s4+s5 + noisamp*randn(size(s1))); %test signal with noise
            dbx1 = dbt(x1,Fs,bandwidth,'lowpass',slow_lowpass,'upsampleFx',2,'remodphase',true,'upsampleTx',bandwidth_ratio-1);
            dbx2 = dbt(x1,Fs,bandwidth*bandwidth_ratio,'lowpass',200,'remodphase',true,'upsampleFx',2);
            dbx = [dbx1,dbx2];

        case {2,'bivariate'}

            %%% Bivariate PAC: x1 contains s1 and x2 s4 + s5
            %%% In this case the fast oscillation envelope in x2 is modulated
            %%% by the slow oscillation in x1. Bivariate bicoherence is
            %%% symmetric under combined sign reversal of the slow frequency
            %%% and permutation of the two fast bands (i.e. f2 and -f2-f1)
            x1 = real(s1+noisamp*randn(size(s2))); %test signal 1 with noise
            x2 = real(s4 + s5+noisamp*randn(size(s2))); %test signal 2 with noise

            dbx1 = dbt(x1,Fs,bandwidth,'lowpass',slow_lowpass,'upsampleFx',2,'remodphase',true,'upsampleTx',bandwidth_ratio-1);
            dbx2 = dbt(x2,Fs,bandwidth*bandwidth_ratio,'lowpass',200,'remodphase',true,'upsampleFx',2);
            dbx = [dbx1,dbx2];


        case {3,'trivariate'}

             %%% Trivariate cross-bicoherence: x1 contains s1, x2 s4, and x3 s5
             %%% Trivariate cross-frequency coupling means that coherency
             %%% between x2 and x3 is modulated by the slow oscillation in
             %%% x1. In this case, the phase difference between x2 and
             %%% x3 follows the phase of the slow oscillation.
             %%% Note that there are two possible sidebands: at fast + slow
             %%% frequency and fast - slow frequency. For this reason, there
             %%% is no necessary symmetry in the trivariate bispectrrum under
             %%% sign reversal or permutation of any frequency axes,
             %%% except in the case of real-valued signals, conjugate symmetry
             %%% under sign reversal of all axes.
             
            x1 = real(s1+noisamp*randn(size(s2))); %test signal 1 with noise
            x2 = real(s4+noisamp*randn(size(s2))); %test signal 2 with noise
            x3 = real(s5+noisamp*randn(size(s2))); %test signal 3 with noise

            dbx1 = dbt(x1,Fs,bandwidth,'lowpass',slow_lowpass,'upsampleFx',2,'remodphase',true,'upsampleTx',bandwidth_ratio-1);
            dbx2 = dbt(x2,Fs,bandwidth*bandwidth_ratio,'lowpass',200,'remodphase',true,'upsampleFx',2);
            dbx3 = dbt(x3,Fs,bandwidth*bandwidth_ratio,'lowpass',200,'remodphase',true,'upsampleFx',2);

            dbx = [dbx1,dbx2,dbx3];
        otherwise
            error('simulate must be ''univariate'', ''bivariate'' or ''trivariate.''')
    end

    %% Compute cross-bispectrum with pspect2
    %%% By default, the bispectrum is computed with respect to f1 and f2,
    %%% where -f1-f2 the implicit third frequeny. With the "symmetrize" option
    %%% it is computed with resepct to f1 and f2-.5*f1, such that .5f1-f2 the implicit 
    %%% 3rd frequency. With this option, the fast oscillation axis in the plot refers to 
    %%% the center of a fast-oscillation band of bandwidth f1, rather than its lower edge.
    %%% The "full_range" option retains negative frequencies in the plot,
    %%% which is used here to illustrate differing symmetries of
    %%% univariate, bivariate and trivariate bicoherence.
    psp = pspect2(dbx,3,'full_range',true,'symmetrize',true,'stats',true);  
    pval = psp.stats.pval;
    pval(end+1) = nan;
    q = fdr(pval);
    PV = pval(psp.squaremat); %Pvalue for test on mean
    Q = q(psp.squaremat);     %FDR q value
    
    subplot(1,length(simulate),k)
    
    BC = abs(psp.pspect)./psp.normalization - psp.bias; %Bias corrected bicoherence 
    imagesc(psp.fs{:},abs(BC)')
    axis xy
    xlabel('Chan. 1: Slow oscillation freq. (Hz)')
        
     switch simulate{k}
        case {1,2,'trivariate','bivariate'}
            title(sprintf('%s cross-bicoherence',simulate{k}))
            if strcmp(simulate{k},'trivariate')
                ylabel('Chans. 2/3: Fast oscillation freq. (Hz)')
                fprintf('\nTrivariate cross-bicoherence is symmetric only with respect to sign reversal of all frequencies.\n')
            else
                ylabel('Chan. 2: Fast oscillation freq. (Hz)') 
                fprintf('\nBivariate cross-bicoherence is symmetric with respect to permutation of fast oscillation frequencies and sign reversal\nof the slow oscillation.\n')

            end
         otherwise
             title(sprintf('%s bicoherence',simulate{k}))
            ylabel('Chan. 1: Fast oscillation freq. (Hz)')
          fprintf('\nUnivariate bicoherence is symmetric with respect to all permutations of frequencies and sign reversal.\n')

     end
     caxis([0 .75])
     drawnow
end
c = colorbar('Position',[0.95    0.2    0.01    0.6]);
ylabel(c,'bicoherence')

