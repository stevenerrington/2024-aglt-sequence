function [B,cv] = blurNorm(X,M,n)

%function B = blurNorm(X,M,n)
%
% Normalizes blur of image X to that of M. 
% BlurNorm Takes the radial average of the log 2D fft of M, then
% subtracts that of X. Ifft is applied to obtain an estimated blurring
% kernel which is then convolved with X to obtain Z. n is the size of the
% the blurring kernel window (default is 10);
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/imagetools/blurNorm.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

if strcmp(class(M),'uint8')
    M = double(M)./255;
end
if strcmp(class(X),'uint8')
    X = double(X)./255;
end

if nargin <3
    n = 9;
end

%Crop M to square.
sz = size(M);
ln = min(sz(1:2));
M = M([1:ln] + floor((sz(1)-ln)./2),[1:ln] + floor((sz(2)-ln)./2),:); 

sz = size(M); sz = sz(1:2);


%%%%%%%%%% Make transformation from linear to radial coordinates
angles = 0:2*pi./sz(1):pi./2;

rs = 0:min(sz./2);

Rs = ones(length(angles),1)*rs;
As = angles'*ones(1,length(rs));

Q1i = Rs.*cos(As)+1;
Q1j = Rs.*sin(As)+1;

Q2i = Rs.*cos(As+pi./2) + sz(1);
Q2j = Rs.*sin(As+pi/2) + 1;
Q3i = Rs.*cos(As+pi) + sz(1);
Q3j = Rs.*sin(As+pi)+ sz(2);
Q4i = Rs.*cos(As-pi./2) + 1;
Q4j = Rs.*sin(As-pi/2) + sz(2);


It = round(cat(1,Q1i,Q2i,Q3i,Q4i));
Jt = round(cat(1,Q1j,Q2j,Q3j,Q4j));

Rota = zeros(sz);
Rota = sub2ind(size(Rota),It,Jt);
%%%%%%%%%%%%%%

 for c = 1:size(M,3)    
        lFm(:,:,c) = log(abs(fft2(M(:,:,c))));
 end
    
lFm = mean(lFm,3);

Y = X;

X = X([1:ln] + floor((sz(1)-ln)./2),[1:ln] + floor((sz(2)-ln)./2),:); 
       
for c = 1:size(X,3)
      lF(:,:,c) = log(abs(fft2(X(:,:,c))));
end
 
lF = mean(lF,3);
RDiff = lFm(Rota)  - lF(Rota);
RDM = mean(RDiff);
RDM = [RDM,fliplr(RDM)];
cv = real(ifft(exp(RDM)));
cv = cv([end-floor(n/2)+1:end,1:ceil(n/2)]);
for c = 1:size(Y,3)
    B(:,:,c) = conv2(cv,cv,Y(:,:,c),'same');
end 
B = B*mean(mean(mean(M)))./mean(mean(mean(B)));
%B(B>1) = 1;
%B(B<0) = 0;

