function Y = lpfilt(X,filt,sband)

%Lowpass filter using finite impulse response function

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin<3
    sband = .2;
end
sband = sband./(filt(1)./2);

if size(X,1) == 1
    X = X';
end


n = size(X,1);

fnorm  = filt(2)./(filt(1)./2);
%mx = mean(X);
%X = X-mx(ones(size(X,1),1),:);
%B = fir1(size(X,1)-1, fnorm);
B = fir2(size(X,1)-1, [0, fnorm(1), fnorm(1)+sband, 1], [1 1 0 0] );

F = abs(fft(cat(2,B,  zeros(1,size(X,1) ))))';

Y = real( ifft( fft(cat(1,X,zeros(size(X)))).*F(:,ones(size(X,2),1)) ) );

Y = Y(1:n,:);

