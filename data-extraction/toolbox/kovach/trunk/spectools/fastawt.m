function [A,frqs,AV,flatten, Compress,F] = fastawt(X,per,freq,sigf,sigt)

%[A,freq,AV,flatten,Compress] = fastawt(X,per,[fs,flo,fhi] ,sigt,sigf);
%Computes Gabor wavelet transform between flo and fhi in logarithmic
%frequency steps sigf and sigt times the time and frequency standard widths, 
%respectively. Per is the center width of the wavelet in units of center period. 
%
% X - signal
% per - period width of wavelet;
% freq - [sf flo fhi]  Sampling freq, and frequency range.
% sigf - step size in frequency  (default 1)
% AV - reduced transform
% flatten - Matrix to reconstruct reduced awt AV(flatten)
% Compress- Sparse matrix allowing compression of full awt,  AV = Compress*A; 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/fastawt.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

X = X(:);

if nargin < 4
    sigf = 1;
end
if nargin < 5
    sigt = 1;
end


sf = freq(1);
nyqf = sf./2;
flo = freq(2)./nyqf; %In units of Nyquist frequency
fhi = freq(3)./nyqf; 
N = length(X);


    %Frequency steps in radian units.
w = exp([ log(flo*2*pi) : sigf./(4*pi*per) : log(fhi*2*pi) ]);


omegas = [0:1/N:1,1/N-1:1/N:-1/N]'.*2*pi;

FX = fft([X;zeros(size(X))]);


% width of wavelets in the frequency domain 
sw = ones(2.*N,1)*w./(2*per*pi); 

% Fourier domain gaussian wavelets. 
F = exp( -((omegas*ones(1,length(w)) - ones(2.*N,1)*w)./(sqrt(2).*sw)).^2);
F = sqrt(2*N)*F./(ones(2.*N,1)*sqrt(sum(F.^2)));


F(N+1:end,:) = 0;

%Convolution of wavelet and signal (via multiplication of Fourier transforms)
a = FX(:,ones(1,length(w))).*F;


A = ifft(a); %inverse fft gives signal convolved with wavelets

A = A(1:N,:)'; %Because the convolution is circular we padded with N points
   

%%%%

if nargout > 2 
    flatten =ones(size(A));
    AV = [];
    ind = 1;
    lng = size(A,2);
    
    if nargout >= 5
           Compress = spalloc(length(AV),prod(size(A)), prod(size(A)));
    end
    
    for i = 1:length(w)

        ds = floor(sigt./(sw(1,i)));

        for t = 1:ds:lng     

            %if strcmp(version,'6.5.0.180913a (R13)')
            %    warning off MATLAB:divideByZero
            %end

            AV(:,ind) = mean(A(i,t:min([t+ds-1,lng])),2);

            flatten(i,t:min([t+ds-1,lng])) = ind;
            
            if nargout == 5
                 Compress(ind,i + size(A,1)*([t:min([t+ds-1,lng])] - 1)) = 1./length( t:min([t+ds-1,lng]) );
            end
            ind = ind+1;    
        end
    end

    AV =squeeze(AV);
    
    

end

frqs = w*sf./(4*pi);