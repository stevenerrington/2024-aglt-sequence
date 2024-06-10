

function [spec,bstr,F1] = mtspect(X1,varargin)


% S = mtspect(X)
%Multitaper spectrum of a signal
%If X is a matrix then the output is computed for each column and
%then averaged across columns.
%
%Options: 
%  nw -  set bandwidth parameter (default = 4)
%
%  ntapers -  set number of tapers (default = 2*mw-1
%   
%  type   -  'dpss' for slepian tapers, 'sine' for sine tapers.
%  
%  bootstrap -  'true' computes confidence intervals with bootstraping.
%               The argument can also be a struct, bstr, with 
%                    bstr.nbstr = number of bootstraps
%                    bstr.ps = confidence levels to compute
%
%See also MTF and MTCSPECT

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/mtspect.m $
% $Revision: 518 $
% $Date: 2014-11-04 18:52:49 -0600 (Tue, 04 Nov 2014) $
% $Author: ckovach $
% ------------------------------------------------

%C Kovach 2008


nw = 4;
ntapers = [];
bootstrap = false;
bstr.nbstr = 250;  %Number of bootstrap samples
bstr.ps = [.005 .01 .05 .25 .5 .75 .95 .99 .995]'; %Confidence levels to return in the bootstrap
type = 'dpss';

i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
       case 'nw'
            nw = varargin{i+1};
            i = i+1;
       case 'ntapers'
            ntapers = varargin{i+1};
            i = i+1;
       case 'type'
            type = varargin{i+1};
            i = i+1;
       case 'bootstrap'
            bootstrap = varargin{i+1};
            if isstruct(bootstrap)
                bstr = bootstrap;
                bootstrap = true;
            end
            
            i = i+1;
            
        otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end
            


nt = size(X1,1);
nc = size(X1,2);

if isempty(ntapers)
    ntapers = 2*nw-1;
end

switch lower(type)
    case {'dpss','slepian'}   %dpss tapers
        
        
        [mtwins, ev] = dpss(nt,nw,ntapers);
%         evs = repmat(ev,nc,1);
    case {'sine','sin'}     % sine tapers
        wx = (1:nt)'./(nt+1)*pi;
        ks = 1:ntapers;
        mtwins = sin(wx*(ks))*sqrt(2/(nt+1));
        ev = ones(ntapers,1);
    otherwise
        
        error('%s is not a recognized family of tapers.',type)
end

evs = repmat(ev,nc,1);
        
        
W = repmat(mtwins,1,nc);
X1 = kron(X1,ones(1,ntapers));
% X2 = kron(X2,ones(1,ntapers));

F1 = fft(W.*X1);
% F2 = fft(W.*X2);


P1 = F1.*conj(F1);
% P2 = F2.*conj(F2);
% P12 = F1.*conj(F2);

% cspec = P12*evs;
spec = P1*evs;
% specs(:,2) = P2*evs;

% coh = abs(cspec)./sqrt(prod(specs,2));


if bootstrap ==1
    
    np1 = fprintf('Bootstrap...');
    bs = '\b';
    bstr.cspecbs = zeros(size(F1,1),bstr.nbstr);
    for j = 1:bstr.nbstr
%         fprintf('%03i',j);
        shuff = kron((ceil(rand(1,nc)*nc)-1)*ntapers, ones(1,ntapers)) + kron(ones(1,nc),1:ntapers);
        
%        cspecbs = P12(:,shuff)*evs;
        
        specsbs(:,1) = P1(:,shuff)*evs;
%         specsbs(:,2) = P2(:,shuff)*evs;
%         bstr.cohbstr(:,j) = abs(cspecbs)./sqrt(prod(specsbs,2));
         bstr.specsbs(:,j) = specsbs;

%         fprintf('\b\b\b');
    end
%         bstr.cspecbs = cspecbs;
        bstr.qtile = quantile(bstr.specsbs',bstr.ps);
        fprintf(repmat(bs,1,np1));

% elseif  bootstrap == 2  %do jacknife
%     
%     [bstr.jkm,bstr.jkse] = cohjacknife(F1,F2,evs);
    
    
  
else
    bstr = [];
end





