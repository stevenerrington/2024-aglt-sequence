
function [xrs,fsrs]  = ftresamp(x,fs,targetfs,targetn,tol)

%xrs = ftresamp(x,fs,targetfs,targetn)
%Resample x to a sampling rate of targetfs and pad/truncate duration
%of targetn using Fourier and symmetric time-domain padding/truncation. 


n = size(x,1);

if nargin < 4 || isempty(targetn)
    targetn = round(length(x)*targetfs/fs);
end

if nargin < 5 || isempty(tol)
    %%% Set error tolerance to one sample or less within the observation window.
    tol = 1/targetn;
end

%%% Pad x 
x(2*n,:) = 0;
x = circshift(x,ceil(n/2));
n = 2*n;


[p,q] = rat(targetfs/fs,tol);
fsrs = fs*p/q;


fx = fftshift(fft(x),1);
if targetfs<fs % Add a lowpass filter
   w = ((0:length(fx)-1)- floor(length(fx)/2))/length(fx)*fs;
   fx(abs(w)>fsrs*(1/2-1/n),:)=fx(abs(w)>fsrs*(1/2-1/n),:)*.5;
   fx(abs(w)>fsrs/2,:)=0;
end

%%% Fourier-domain upsampling by a factor of p
if p > 1
    fx(p*n,:) = 0;
end
fx = circshift(fx,-floor(n/2));

xrs = real(ifft(fx))*p;

%%% Add an offset that aligns the centermost samples (so that interpolation error is minimized there)
offs = mod(floor(length(xrs)/2),q);

xrs = xrs(1+offs:q:end,:); %Decimate by q

if ~isempty(targetn) %Pad or truncate to specified window length, if given.
    dn = length(xrs)-targetn;
    if dn>0
        xrs = xrs(ceil(dn/2+.5):end-ceil(dn/2),:);
    else
        xrs(targetn,:) = 0;
        xrs = circshift(xrs,floor(-dn/2));
    end
end


