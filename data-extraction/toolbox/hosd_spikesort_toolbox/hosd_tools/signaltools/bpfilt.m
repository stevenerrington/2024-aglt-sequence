function [Y,F] = bpfilt(X,filt,sband,filtord)

%Bandpass filter using finite impulse response function
if nargin < 3 || isempty(sband)
% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

    sband = .2;
end

if size(X,1) == 1
    X = X';
end


sbands = sband.*[1 1]./(filt(1)./2);

if sband(1) >filt(2);
    sbands(1) = filt(2)./(filt(1));
end

if nargin < 4
    filtord = size(X,1)-1;
end

np2 = 2^nextpow2(size(X,1)*2);
X2 = zeros(np2,size(X,2));
n = size(X,1);
X2(1:n,:) = X;

fnorm  = filt(2:3)./(filt(1)./2);
%mx = mean(X);
%X = X-mx(ones(size(X,1),1),:);
%B = fir1(size(X,1)-1, fnorm);
B = fir2(filtord , [0, fnorm(1) - sbands(1), fnorm(1), fnorm(2), fnorm(2) + sbands(2), 1], [0 0 1 1 0 0] );

F = abs(fft(cat(2,B,  zeros(1,size(X2,1)-length(B) ))))';

Y = real( ifft( fft(cat(1,X2)).*F(:,ones(size(X2,2),1)) ) );

Y = Y(1:n,:);

