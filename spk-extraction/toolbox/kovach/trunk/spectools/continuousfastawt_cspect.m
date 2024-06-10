
function [coh,cspect,frqs,raw1,raw2] = continuousfastawt_cspect(X1,X2,freq,per,Tchop)

if nargin < 5
    Tchop = 1:size(X1);
end

fs = freq(1);
[A,frqs,AV,flatten, Compress,F] = fastawt(1/fs:1/fs:10./freq(2)*per,per,freq);

irf = fftshift(ifft(F),1);

fprintf('\n%i banks, doing %3.0i',size(irf,2),0);

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/continuousfastawt_cspect.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

for i = 1:size(irf,2);
    
    fprintf('\b\b\b%3.0i',size(irf,2),i);
    
    cv =zeros(size(X1));
    
    cv(1:size(irf,1)) = irf(:,i);
    
    cv = circshift(cv,-round(size(irf,1)/2));
    
    filt1 = ifft(fft(X1).*conj(fft(cv)));
    filt2 = ifft(fft(X2).*conj(fft(cv)));
    
    cspect(i,:) = mean(filt1(Tchop).*conj(filt2(Tchop)),2);
    coh(i,:) = cspect(i,:)'./sqrt(mean(abs(filt1(Tchop)).^2,2).*mean(abs(filt2(Tchop)).^2,2));
    
    if nargout > 3
        raw1(i,:,:) = filt1(Tchop);
        raw2(i,:,:) = filt2(Tchop);
    end
    
    
end

