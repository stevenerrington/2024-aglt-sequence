function [X,disc,tind,tsort,texpected] = multitrim(X,varargin)

% [X,disc,tind,tsort] = multitrim(X,trim,exvar)
%Removes outliers from X using Multivariate trimming of extreme data points. See Rencher p. 29.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/multitrim.m $
% $Revision: 8 $
% $Date: 2011-05-02 15:39:11 -0500 (Mon, 02 May 2011) $
% $Author: ckovach $
% ------------------------------------------------


trim = .2;
exvar = .9;
pckeep = [];
i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
       case 'trim'
            trim = varargin{i+1};
            i = i+1;
       case 'exvar'
            exvar = varargin{i+1};
            i = i+1;
       case 'pckeep'
            pckeep = varargin{i+1};
            i = i+1;
        otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end


X = squeeze(X);



Xc = X;


mx = mean(Xc,2);
n = size(Xc,2);
%MX = mx(:,ones(size(Y,2),1));

Xc = Xc-mx*ones(1,size(Xc,2));


if size(Xc,1) > size(Xc,2)
    flip = true;
else
    flip = false;
end

if flip
    
    S = (Xc'*Xc)/(n-1);   %Because S is singular, reduce dimensionality with svd.
else
    S = (Xc*Xc')/(n-1);   
end    


totpow = sum(X(:).^2);
if isempty(pckeep)
    D =0
    geteigs = 0;
    while trace(D)<exvar*totpow
        geteigs = 1+geteigs;
        [U,D] = eigs(S,geteigs);
    end
    pckeep = sum(diag(D) <= exvar*totpow);
    [U,D] = eigs(S,pckeep);
    
else
        [U,D] = eigs(S,pckeep);
end        
    
d = diag(D);

if flip
    U = Xc*U*diag((d).^-.5);        %Converts eigenvectors of X'X to nonzero eigenvectors of XX'
end

% if isempty(pckeep)
%     pckeep = min(find(cumsum(d/sum(d))>=exvar));
% end

r = rank(S);

if r == 1
    error('This script can''t handle univariate data.')
end

% if pckeep > r | isempty(pckeep);
%     pckeep = r;
% end
% 
% if pckeep == 1
%     pckeep = 2
% end
    
    
% d = d(1:pckeep);
% U = U(:,1:pckeep);

Xcu = U'*Xc;

zold = inf;
Disc = [];

Xu = Xcu;  
while 1

 Su = Xcu*Xcu'/(size(Xcu,2)-1);
 R = diag(diag(Su).^-.5)*Su*diag(diag(Su).^-.5);
 rs = R(find(triu(R+1e-15,1)));
 
 znew = .5*log((1+rs)./(1-rs));  % Fisher's z transform
  
 if mean(abs(znew-zold)) < .001 %| exvar == 1
     break
 end
 
 zold = znew;
% T = diag(Xcu'*Su^-1*Xcu); 
  T = diag(Xu'*Su^-1*Xu); 

 [tsort,tind] = sort(T);
 disc = tind(end - ceil(trim*length(tind)) + 1:end);
 %Disc = cat(1,Disc,disc);
 %Xcu(:,disc) = [];
 Xcu = Xu; 
 Xcu(:,disc) = [];
 Xcu = Xcu - mean(Xcu,2)*ones(1,size(Xcu,2));
 if trim == 1
     break
 end
end

if ~isempty(Xcu)
  N = size(Xcu,2);
else
  N = size(X,2); 
end

%F = (nu - pckeep +1)/(nu*pckeep)*tsort;
texpected = sort((N - 1)*pckeep/(N - pckeep)*frnd(pckeep,N-pckeep,1,length(tsort)));
   


X(:,disc) = [];
   
    
    


