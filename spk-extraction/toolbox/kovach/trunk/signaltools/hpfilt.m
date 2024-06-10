function Y = hpfilt(X,filt,sband)

%Bandpass filter using finite impulse response function
% sband = .5;
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/signaltools/hpfilt.m $
% $Revision: 583 $
% $Date: 2015-03-31 16:04:39 -0500 (Tue, 31 Mar 2015) $
% $Author: ckovach $
% ------------------------------------------------


if size(X,1) == 1
    X = X';
end

if nargin < 3
        sband = .5./(filt(1)./2);
else
    sband = sband./(filt(1)./2);
end

n = size(X,1);

fnorm  = filt(2)./(filt(1)./2);
%mx = mean(X);
%X = X-mx(ones(size(X,1),1),:);
%B = fir1(size(X,1)-1, fnorm);
B = fir2(size(X,1)-1, [0, fnorm(1) - sband, fnorm(1), 1], [0 0 1 1] );

F = abs(fft(cat(2,B,  zeros(1,size(X,1) ))))';

Y = real( ifft( fft(cat(1,X,zeros(size(B,2),size(X,2)))).*F(:,ones(size(X,2),1)) ) );

Y = Y(1:n,:);

