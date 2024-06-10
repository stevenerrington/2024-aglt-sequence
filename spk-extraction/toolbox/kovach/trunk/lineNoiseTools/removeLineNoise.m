
function [x,noise] = removeLineNoise(x,fs,varargin)

%  [x_denoised,noise] = removeLineNoise(x,fs)
% This script attempts to remove line noise and harmonics by (1)
% Identifying peaks in the spectrum and classifying them as fundamentals or
% harmonics. (2) For each fundamental (F0), computing the FFT of the signal in
% overlapping windows, where window size is as near an integer multiple of
% the fundamental period as possible. (3) For each window, finding the set of harmonics within +/-
% .5 of F0 which have the highest total power and (4) removing them. (5)
% Applying the inverse Fourier transform. (6) Constructing the denoised
% signal from overlapping denoised windows as a weighted average, where weighting as
% a smoothly varying window function. 
%
% See Also REMOVELINENOISE and DENOISE60HZ

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/lineNoiseTools/removeLineNoise.m $
% $Revision: 389 $
% $Date: 2013-10-18 15:53:08 -0500 (Fri, 18 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2010




lothresh = 40; % Only apply denoising to frequencies above this value(!)
fitmin = 20; %minimum frequency used in polynomial fitting

% hithresh = fs/2;
% hithresh = 3000;  %By default Waveclust implements a 3K high pass filter, so we can ignore everything abov 3.
hithresh = min([3000 fs/2]);  %By default Waveclust implements a 3K high pass filter, so we can ignore everything abov 3.
                              % & ignore above nyquist.
% lothresh = 0; % Only apply denoising to frequencies above this value(!)
% hithresh = Inf;  %By default Waveclust implements a 3K high pass filter, so we can ignore everything abov 3.

% if nargin < 3 || isempty(sdthresh)
    sdthresh = 5;  %Extract peaks above this SD 
% end
makeplot = true;
windowSizeRange = [1 4]; % range of acceptable window sizes  
smooth_ord = 8; %Order of smoothing polynomial
repeat = 1; %Repeat denoising this many times times.

use_parallel = false;
% frequencyResolution = .1; %Estimate periodogram with this resolution.


smwin = .2; % Size of the smoothing window applied to the amplitdue of the raw FFT of the signal to estimate spectral power.

i = 1;
while i < length(varargin)
   switch lower(varargin{i})
       case 'freqrange'
           fqr = varargin{i+1};
           lothresh = fqr(1);
           fitmin = fqr(1);
           hithresh = fqr(2);
           i = i+1;
       case 'lothresh'
           lothresh = varargin{i+1};
           i = i+1;
       case 'hithresh'
           hithresh = varargin{i+1};
           i = i+1;
       case 'fitmin'
           fitmin = varargin{i+1};
           i = i+1;
       case 'makeplot'
           makeplot= varargin{i+1};
           i = i+1;
       case 'smoothord'
           smooth_ord = varargin{i+1};
           i = i+1;
       case 'sdthresh'
           sdthresh= varargin{i+1};
           i = i+1;
       case 'windwsizerange'
           windowSizeRange= varargin{i+1};
           i = i+1;
       case 'repeat'
           repeat = varargin{i+1};
           i = i+1;
       otherwise
           error('%s is an unrecognized option.',varargin{i});
   end
   i = i+1;
end
N = length(x);
linecols = 'bcmy';
for rep = 1:repeat
    F = fft(x);
    sm = round((N./fs)*smwin);
    % F = F(1:ceil(N./2));
    P = abs(F);
    g =gausswin(sm*3);
    g = g./sum(g);
    % B = fir2(size(x,1)-1, [0, 8*pi./(smwin*fs)   , 8*pi./(smwin*fs) + 2/fs*10, 1], [1 1 0 0] );
    P = fastconv(P,g); %Smooth P. Alternately we could window in the time domain-- faster.


    % [P,w] = pwelch(x,tukeywin(round(fs/frequencyResolution)),0); %Estimate spectrum


    frq = (0:N-1)./N*fs;

    getfitfrq = frq>fitmin & frq < hithresh;
    fitfrq = frq(getfitfrq);

    getf = frq>lothresh & frq<hithresh;
    frq = frq(getf);


    % % P = convn(P,g,'same');


    mmrescale = @(x) (x-min(x))./(max(x)-min(x));  %A function to Rescale input polynomial estimation to improve numerical stability

    %%% Fit a polynomial to log(P) vs. log(frq) in order to estimate the baseline spectrum
    [ll,~] = polyfit(mmrescale(log(fitfrq)),log(P(getfitfrq))',smooth_ord );
    est = polyval(ll,mmrescale(log(fitfrq)));

    P = P(getf); %Extract desired frequency range
    est = est(getf(getfitfrq));

    %%% Residual after subtracting the baseline
    %%% Roughly, the log ratio of power relative to the surrounding baseline
    resid = log(P)- est'; 

    nanzscore = @(x) (x-nanmean(x))./nanstd(x); % a function to compute the zscore of a vector with nans
    z = nanzscore(resid);
    % sdthresh = norminv(1-sm/length(z));
    ztemp = z;
    rtemp = resid;

    %%% Iteratively lower the threshold, excluding points above threshold until
    %%% none are above
    nabove = sum(abs(z)>sdthresh);
    while nabove > 0 
        rtemp(abs(ztemp)>sdthresh) = nan;
        ztemp = nanzscore(rtemp);
        nabove = sum(abs(ztemp)>sdthresh);
    %     nabove
    end

    %%% Compute z based on adjusted mean and std dev
    z = (resid - nanmean(rtemp))./nanstd(rtemp);


    zthresh = z>sdthresh; % Get values above the specified threshold

    zthresh(end) = 0;
    dzthresh = diff([0;zthresh]);
    pkmin =  find(dzthresh>0); 
    pkmax = find(dzthresh<0);
    pkwidth = pkmax-pkmin; 

    %%% The peak is defined as
    %%% the maximum within the interval from pkmin to pkmax.

    if isempty(pkwidth)
        fprintf('\nNo peaks were above the threshold. No denoising was applied.\n')
        return
    end
     
    linepks = zeros(size(pkwidth));
    for i = 1:length(pkwidth)
       intvl = pkmin(i):pkmax(i);
       [~,mxi] = max(z(intvl));  
       linepks(i) = intvl(mxi);
    end


    % [pk,pks] = getpeak(P);
    % pk = find(pks==1);
    % linepk = pk(z(pk)>sdthresh);
    % % P = fastconv(z,g);


    linefrq = frq(linepks);

    htol = .001; %tolerance for deciding whether a value is an integer multiple of another
    
    harmonics = false(length(linepks));
    for i = 1:length(linefrq)
        htest = linefrq./linefrq(i);
        fpart = htest-round(htest);
        ipart = round(htest);

    %     harmonics(i,:) = abs(fpart) < htol;% & htest~=1;
        harmonics(i,:) = abs(fpart) < htol; % All integer multiples of a given frequency are defined as harmonics


        %%% If multiple closely spaced peaks are around the same frequncy, keep
        %%% only the greater
        samepk = ipart == 1 & harmonics(i,:); 
        pkw = zeros(1,size(harmonics,2));
        pkw(samepk) = pkwidth(samepk);
        [~,mxi]= max(pkw);        
        harmonics(i,mxi) = 0;

    end
    keeppks = ~any(harmonics); %Keep only peaks that are not harmonics  of any other 
    f0s = linefrq(keeppks);    % Fundamental frequencies


    if makeplot
        plot(frq,z,linecols(mod(rep-1,length(linecols))+1))
        hold on,
        plot(frq([1 end]),[1 1]*sdthresh,'g')
        plot(frq(linepks),z(linepks),'r.')
        plot(frq(linepks(keeppks)),z(linepks(keeppks)),'ro')
        axis tight
        title('Extracted Peaks')
        xlabel Hz
        legend({'spectrum','threshold','identified peaks','identified fundamentals'})
        drawnow
    end

    %     wsize = windowsize*fs; 

    %%% Iterate through extracted fundamental frequencies, applying
    %%% denoiseContinuous to each. SEE DENOISECONTINUOUS
    if ~use_parallel
        for i = 1:length(f0s);
            fprintf('\n%i of %i, Removing %0.1f hz line noise',i,length(f0s),f0s(i))
        %     wm = round(wsize*f0s(i)./fs);
            xdn = denoiseContinuous(x,fs,f0s(i)+[-.5 .5],windowSizeRange);
            if nargout > 1
                noise(:,i,rep) = x-xdn;
            end
            x = xdn;
        end
    else
         parfor i = 1:length(f0s);
            fprintf('\n%i of %i, Removing %0.1f hz line noise',i,length(f0s),f0s(i))
        %     wm = round(wsize*f0s(i)./fs);
            xdn = denoiseContinuous(x,fs,f0s(i)+[-.5 .5],windowSizeRange);
            noise(:,i) = x-xdn;
         end
         x = x-sum(noise,2);
    end
end

%%%%%%%%%s

% 
% function [pk,pksign,zeroc] = getpeak(x)
% 
% 
% dx = diff(x);
% sdx = sign(dx);
% 
% dsdx = [0;sdx] - [sdx;0];
% 
% zeroc = sign(x(1:end-1))~=sign(x(2:end));
% % dsdx = diff(sdx);
% 
% pk = abs(dsdx)==2;
% pk([1 end]) = true;
% 
% pksign = sign(dsdx.*pk);


%%%%%%%%%%%%

function cab = fastconv(A,B)

nb = length(B);

% front = 1:ceil(nb/2);
% back = ceil(nb/2)+1:nb;

ab = zeros(size(A));

ab(1:ceil(nb/2)) = B(floor(nb/2)+1:end);
ab(end - floor(nb/2)+1:end) = B(1:floor(nb/2));

cab = ifft(fft(A).*conj(fft(ab)));
