function [A,frqs,AV,flatten, Compress,F] = continuousfastawt(X,per,freq,sigf,sigt)

%[A,freq,AV,flatten,Compress] = fastawt(X,per,[fs,flo,fhi] ,sigt,sigf);
%Computes Gabor wavelet transform between flo and fhi in logarithmic
%frequency steps sigf and sigt times the time and frequency standard widths, 
%respectively. Per is the center width of the wavelet in units of center period. 
%This version is for continuous signals and uses no padding to improve efficiency.
%
% X - signal
% per - period width of wavelet;
% freq - [sf flo fhi]  Sampling freq, and frequency range.
% sigf - step size in frequency  (default 1)
% AV - reduced transform
% flatten - Matrix to reconstruct reduced awt AV(flatten)
% Compress- Sparse matrix allowing compression of full awt,  AV = Compress*A; 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/continuousfastawt.m $
% $Revision: 115 $
% $Date: 2012-05-24 23:50:07 -0500 (Thu, 24 May 2012) $
% $Author: ckovach $
% ------------------------------------------------

% Written by C. Kovach




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
omegas = omegas(1:2:end);

FX = fft(X);

% width of wavelets in the frequency domain 
sw = repmat(w,N,1)./(2*per*pi); 

% Fourier domain gaussian wavelets. 
F = exp( -((omegas*ones(1,length(w)) - ones(N,1)*w)./(sqrt(2).*sw)).^2);
F = sqrt(N)*F./(ones(N,1)*sqrt(sum(F.^2)));


F(N+1:end,:) = 0;  %Zero out the upper frequencies to get a complex signal
F(2:N-1,:) = 2*F(2:N-1,:);  % This is equivalent to applying hilbert transform to filtered data
if mod(N,2) == 0
    F(N,:) = 2*F(N,:);
end
    

% Convolution of wavelet and signal (via multiplication of Fourier transforms)
% For efficiency, no padding is used. There will be leakage at the edges, but 
% we are assuming that the signal is a long continuous recoring for which  edge effects are  
% a relatively minor problem.

a = FX(:,ones(1,length(w))).*F;


A = ifft(a); %inverse fft gives signal convolved with wavelets

% A = A(1:N,:)'; %Because the convolution is circular we padded with N points
   

%%%%

if nargout > 2
    
    flatten =ones(size(A));
    AV = [];
    ind = 1;
    lng = size(A,2);
    
    if nargout >= 5
           Compress = spalloc(length(AV),numel(A), numel(A));
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