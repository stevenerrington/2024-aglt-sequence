
function xcr = xcorrnorm(aa,bb)

%function xcr = xcorrnorm(X,Y)
%
% Sliding normalized cross correlation cross correlation. This gives the 
% correlation between column vector X and a sliding window in column Y.
% It differs from a simple cross correlation in that the normalization
% factor at each time point includes only the window of Y which overlaps
% X.
% 
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C. Kovach 2012

npts = max(length(aa),length(bb))*2;
aapad = detrend(aa,'constant');
aapad(end+1:npts) = 0;
bbpad = detrend(bb,'constant');
bbpad(end+1:npts) = 0;

onespad = ones(size(aa));
onespad(end+1:npts) = 0;

bmn = ifft(conj(fft(onespad)).*fft(bbpad))/sum(onespad);
amn = ifft(conj(fft(onespad)).*fft(aapad))/sum(onespad);
xcv = ifft(conj(fft(aapad-amn)).*fft(bbpad-bmn));
asq = sum(aapad.^2);
bsq = ifft(conj(fft(onespad)).*fft((bbpad-bmn).^2));
xcr = xcv./sqrt(abs(bsq.*asq )+eps);
% xcr = xcr(1:length(bb),:);