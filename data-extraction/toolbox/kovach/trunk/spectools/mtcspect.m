

function [coh,cspec,specs,bstr,F1,F2] = mtcspect(X1,X2, nw,bootstrap,bstr)

%Multitaper coherence and crosspectrum for inputs X1 and X2
%If X1 and X2 are matrices then the output is computed for each column and
%then averaged across columns(!).
%
%See also MTF

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/mtcspect.m $
% $Revision: 177 $
% $Date: 2013-03-13 18:00:59 -0500 (Wed, 13 Mar 2013) $
% $Author: ckovach $
% ------------------------------------------------

if nargin < 4
    bootstrap = false;
end
if nargin <5 && ismember(bootstrap,[1 3])
    bstr.nbstr = 250;
    bstr.ps = [.005 .01 .05 .25 .5 .75 .95 .99 .995]';
end


%C Kovach 2008
nt = size(X1,1);
nc = size(X1,2);

ntapers = 2*nw-1;
[mtwins, ev] = dpss(nt,nw,ntapers);

W = repmat(mtwins,1,nc);
X1 = kron(X1,ones(1,ntapers));
X2 = kron(X2,ones(1,ntapers));
evs = repmat(ev,nc,1);

F1 = fft(W.*X1);
F2 = fft(W.*X2);


P1 = F1.*conj(F1);
P2 = F2.*conj(F2);
P12 = F1.*conj(F2);

cspec = P12*evs;
specs(:,1) = P1*evs;
specs(:,2) = P2*evs;

coh = abs(cspec)./sqrt(prod(specs,2));


switch bootstrap
    case 1
    
    np1 = fprintf('Bootstrap...');
    bs = '\b';
    bstr.cspecbs = zeros(size(F1,1),bstr.nbstr);
    for j = 1:bstr.nbstr
%         fprintf('%03i',j);
        shuff = kron((ceil(rand(1,nc)*nc)-1)*ntapers, ones(1,ntapers)) + kron(ones(1,nc),1:ntapers);
        
       cspecbs = P12(:,shuff)*evs;
        
        specsbs(:,1) = P1(:,shuff)*evs;
        specsbs(:,2) = P2(:,shuff)*evs;
        bstr.cohbstr(:,j) = abs(cspecbs)./sqrt(prod(specsbs,2));

%         fprintf('\b\b\b');
    end
%         bstr.cspecbs = cspecbs;
        bstr.qtile = quantile(bstr.cohbstr',bstr.ps);
        fprintf(repmat(bs,1,np1));

    case 2  %do jacknife
        fprintf('\nJackknife...')
        [bstr.jkm,bstr.jkse] = cohjacknife(F1,F2,evs);
    
    
    case 3 % permutation test
        np1 = fprintf('Permutation test...');
        bs = '\b';
        bstr.cspecbs = zeros(size(F1,1),bstr.nbstr);
        specsbs(:,2) = P2*evs;

        for j = 1:bstr.nbstr
    %         fprintf('%03i',j);
            shuff = kron((randperm(nc)-1)*ntapers, ones(1,ntapers)) + kron(ones(1,nc),1:ntapers);
            
           P12 = F1(:,shuff).*conj(F2);
           cspecbs = P12*evs;

            specsbs(:,1) = P1(:,shuff)*evs;
%             specsbs(:,2) = P2(:,shuff)*evs;
            bstr.cohbstr(:,j) = abs(cspecbs)./sqrt(prod(specsbs,2));

    %         fprintf('\b\b\b');
        end
%         bstr.cspecbs = cspecbs;
        bstr.qtile = quantile(bstr.cohbstr',bstr.ps);
        fprintf(repmat(bs,1,np1));

        
    otherwise
    bstr = [];
end





