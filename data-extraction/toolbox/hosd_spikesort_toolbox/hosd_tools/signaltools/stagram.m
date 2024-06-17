function [S,Nsp,t,wt,SW] = stagram(spikes,signal,windowsize,stepsize,avgwin,Fs)

%  [S,Nsp,t,wt] = stagram(spikes,signal,windowsize,stepsize,avgwin,Fs)
%
% function stagram(spikes,signal,windowsize)
%
% Computes a moving spike triggered sum over the window of the given
% size.
%
% The columns of S are the STA for spikes which occur within windows of size 
% windowsize centered at t. wt is the time index for each row of S.
% 
% avgwin is the size of the window in seconds (or samples for Fs=1) over
% which to average (determining the number of rows in S). It can be a
% single value w, for which the range is [-w/2,w/2], or a two element vector
% giving a time range relative to the spike. Note that avgwin must include
% the window on the spike train and should include padding at least as
% great as  windowsize/2.
%
% Nsp is the number of spikes observed in each window.
%


% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C Kovach 2012

normalize_window_bias = true; %%% Convolve the signal with the window function 
                                %%% and subtract to account for bias.

if nargin < 3 || isempty(windowsize)
    windowsize=256;
end
if nargin < 4 || isempty(stepsize)
    stepsize = windowsize/2;
end
if nargin < 6 || isempty(Fs)
    Fs = 1;
end

Nw = ceil(length(signal)/windowsize);
Nt = ceil(length(signal)/windowsize)*windowsize;
Ntw = Nt/Nw;
Nstep = Nt/stepsize;

if nargin < 6 || isempty(avgwin)
    avgwin = Ntw./Fs;
end
if length(avgwin) == 1
   avgwin = avgwin*[-1 1]; 
end
    
Nsw=windowsize/stepsize;

if Nt > length(signal)
    signal(Nt)=0;
end
if Nt > length(spikes)
    spikes(Nt)=0;
end


Navwin = round(avgwin*Fs);
awin = Navwin(1):Navwin(2);
Naw = length(awin);

%%% reshape

ind = repmat(awin',1,Nstep-2*Nsw) + repmat((Nsw+1:Nstep-Nsw)*stepsize,Naw,1);
ind=ind(:,~any(ind < 1 | ind > Nt));

spikes=spikes(ind);

W = 0*spikes;
W(awin >= -Ntw/2 & awin <= Ntw/2,:)=1;
%%% Pad spikes so that the spike average extends outside of the window
spikes=spikes.*W;

signal = signal(ind);

% if window_signal 
%     signal(awin < -Ntw/2 | awin > Ntw/2,:)=0;
% end

%%% FFT convolution

Fsig = fft(signal);
Fspike = fft(spikes);


%%% Take average
Nsp = sum(spikes);

if normalize_window_bias 
    FW = fft(W).*repmat(Nsp./sum(W)+eps,Naw,1);
%     FW = fft(W);%.*repmat(Nsp./sum(W)+eps,Naw,1);
    SW= ifft(Fsig.*conj(FW));
    S = ifft(Fsig.*conj(Fspike-FW));
%     if normalize_window_bias;
%         S=S./;
%     end
else
    S = ifft(Fsig.*conj(Fspike));
    SW=0;
end
S = fftshift(S,1)./ repmat(Nsp+eps,Naw,1);
if nargout > 4
    SW = fftshift(SW,1);
end
% S = fftshift(ifft( fft(signal).*conj(fft(spikes))),1);
% WS = fftshift(ifft( fft(signal).*conj(fft(W))),1);



t = ind(awin==0,:)./Fs;
wt = awin./Fs;



 



