function denoised = denoise60hz(x,fs,frange)

% Line noise removal based on identification of fundamental and harmonic components of a periodic signal.
%
% x - a matrix of data with rows representing time and columns representing trials
% fs - sampling frequency
%  frange - range of frequencies to search for line noise. If not specified,
%           by default it is [59.5 60.5]. Within this range of fundamental frequencies 
%           the value the value which maximizes the total power of fundamental 
%           and harmonic fourier components is extracted. 
%
% NB, To work properly, the duration of x should be an integer multiple
% of the fundamental frequency of the line noise.
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/lineNoiseTools/denoise60Hz.m $
% $Revision: 128 $
% $Date: 2012-06-08 20:37:23 -0500 (Fri, 08 Jun 2012) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2008

if nargin < 3
    frange = [59.5 60.5]; %Range of frequencies over which to look
end
step = .01; %Size of frequency step within the range
N = size(x,1);




nfrange = frange./(fs./2); %Converting frequencies to normalized 
nstep = step./(fs./2);


testfreq = nfrange(1):nstep:nfrange(2); % Looks at power in fundamental and harmonics of these frequencies

cwidths = testfreq.*N./2; %Frequency steps in terms of sampling units

%The following constructs a set of dirac combs which correspond to 
% periodic signals with fundamental frequencies in testfreq

Combs = zeros(size(x,1),length(cwidths));

for i = 1:length(cwidths)
    
    Combs(round([cwidths(i):cwidths(i):N./2]+1),i) = 1;
    Combs(round(end - [cwidths(i):cwidths(i):N./2] + 1),i) = 1;
        
end

F = fft(x);

pow = abs(F).^2'*Combs;  %Amount of power explained by each comb 

[mx,maxcomb] = max(pow'); %Locating maximum for each trial

noiseft = F.*Combs(:,maxcomb); %Noise is extracted by multiplying the fft of the signal by the best comb
                                   
noise = ifft(noiseft);

denoised = x - noise; %denoised signal
