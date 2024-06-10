function [parest,err,lik,S] = tlsmle(X)

% [parest,err,lik,S] = tlsmle(X)
% Finds the maximum likelihood parameter estimates for a location-scale
% t-distribution using Newton's method.
%
% Returns a vector parest = [mu sigma nu], where mu is the location parameter,
% sigma the scale paremeter and nu is the shape parameter.
%
% If X is a matrix, data are fitted separately for each column of X.
%
% Err is asymptotic estimation error and S is the asymptotic error covariance.

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/tlsmle.m $
% $Revision: 504 $
% $Date: 2014-05-25 00:11:58 -0500 (Sun, 25 May 2014) $
% $Author: ckovach $
% ------------------------------------------------

% C Kovach 2014

%%% starting estimates
parest(1,:) = median(X); % mu
parest(3,:) = 6./kurtosis(X)+4; %nu
parest(2,:) = 1./(parest(3,:).*var(X)); % For the purposes of fitting this will be parameterized as h/nu
                                     % Starting here to ensure that 2nd
                                     % derivative is convex.
%parest(3,:) = 3;  % nu

tol = 1e-6;

derr = Inf;

N = size(X,1);
ncol = size(X,2);

%%% error term
err = @(par)(X-repmat(par(1,:),N,1));  

%%% error squared times precision/shape
Zsq = @(par)  err(par).^2*diag(par(2,:)); 

%%%Difference of digamma functions, which appears in the likelihood derivatives 
dpsi = @(par,k)  N/2^(k+1)*(psi(k,(par(3,:)+1)/2)-psi(k,par(3,:)/2)); 

%%% Likelihood function
L = @(par) N* ( .5*( log(par(2,:))  - log(pi) ) +...
          gammaln((par(3,:)+1)/2) - gammaln(par(3,:)/2) )...
           - (par(3,:)+1)/2.*sum(log(1+Zsq(par)));

%%% Likelihood derivative
DLDmu = @(par)  par(2,:).*(par(3,:)+1).*sum(err(par)./(1+Zsq(par)));
DLDh = @(par)  .5*N./par(2,:) - .5*(par(3,:)+1).*sum(err(par).^2./(1+Zsq(par)));
DLDnu = @(par) dpsi(par,0) - .5*sum(log(1+Zsq(par))) ;


%%% Second derivatives
D2LDmu2 = @(par) (par(3,:)+1).*par(2,:).*sum((Zsq(par)-1)./(1+Zsq(par)));                
D2LDh2 = @(par)  (par(3,:)+1)/2.*sum(err(par).^4./(1+Zsq(par)).^2) -.5*N./par(2,:).^2 ;                
D2LDnu2 = @(par) dpsi(par,1);
D2LDhDmu = @(par) par(2,:).*(par(3,:)+1).*sum(err(par)./(1+Zsq(par)).^2);
D2LDnuDmu = @(par) par(2,:).*sum(err(par)./(1+Zsq(par)));
D2LDnuDh = @(par)  - .5*sum(err(par).^2./(1+Zsq(par)));


DL = @(par) sparseblock(cat(1,DLDmu(par),DLDh(par),DLDnu(par)),1);

D2L = @(par) sparseblock(reshape(cat(1,D2LDmu2(par),D2LDhDmu(par),D2LDnuDmu(par),D2LDhDmu(par),D2LDh2(par),D2LDnuDh(par),D2LDnuDmu(par),D2LDnuDh(par),D2LDnu2(par)),3,ncol*3),3);


iter = 0;
maxiter = 1000;
maxnu = 500;
lvmult = 1;
lvm = exp(lvmult*4).*ones(3,ncol);
ll = L(parest);
dstep = zeros(size(parest));
d2step = zeros(size(parest));
regH = 0;

while derr > tol && iter < maxiter
    
    
    H = D2L(parest) ;
   % iD = -iD.*diag(sign(diag(iD)));
    dstep(:) = -(H-diag(lvm(:)))\unsparsify(DL(parest),'transpose');
    
    %%% Second derivative in the direction of the step should be negative 
    %%% if we are converging on a maximum, hence flip the sign if it is
    %%% positive.
    crv = sum(H*dstep.*dstep);
    if any(crv>0) 
        modH = [1 1 0 ]'*(crv>0);
        H = H - 2*(modH(:)*modH(:)').*H;
        dstep(:) = -H\unsparsify(DL(parest),'transpose');
 
        0;
    end
    
   
    %parestold = parest;
    parest = parest + dstep;

    if any(parest(3,:) > maxnu)
        clp = diag([0 0 1])*parest>maxnu;
        parest(clp)=maxnu;
        dstep(clp)=0;       
    end
      derr = sqrt(sum(dstep.^2)./ncol);
  
     llnew = L(parest);
     ww = llnew <= ll;
     if any(ww)
         lvm = exp(lvm + lvmult*[1 1 1]'*ww);
        parest = parest - dstep*diag(ww);
     end
     ll(~ww) = llnew(~ww);
      lvm = exp(log(lvm) - lvmult*[1 1 1]'*(~ww));
         
         
         
%     llnew = L(parest);
%     wreg = (crv > 0||llnew<=ll);
%     if any(wreg)
%         
%         ww = sparseblock(dstep,1);
%         
%         regH = dstep(:)'*diag(exp(lvm(:) + wreg(:)*lvmult))*dstep(:);
        
%     if any(ll < llold)
%         
%         ddir = dstep(:)'*DL(parest);
%         d2step(:) =D2L(parest)*dstep(:);
%         d2dir = sqrt(sum(dstep.*d2step)./lvmmult); %% Second derivative in the direction of the step (should be negative if we are converging on a maximum)
%         fix = d2dir > 0;
%         spds =  sparseblock(dstep*diag(d2dir.*fix),1);
%         
%         
%         lvm = lvm - 2*(spds'*spds); %Flip the curvature if 2nd derivative is positive.
%         
%         parest(:,fix) = parestold(:,fix);
%         
%     else
%         lvm = 
        
    
        
    iter = iter+1;
    
    lls(iter) = L(parest);
end

if iter >= maxiter
    warning('Failed to converge after %i iterations.\nLast step size = %0.2g',maxiter,full(derr))
end
if any(parest(3,:)==maxnu)
    warning('Shape parameter reached an upper bound of %i. These data are probably not leptokurtic.',maxnu)
end

if nargout > 1
    %Contrast to obtain Hessian for the reparameterized sigma 
    cm = [ones(size(parest(1,:))); zeros(size(parest(1,:))); zeros(size(parest(1,:)))];
    cs = [zeros(size(parest(1,:))); -.5./sqrt(parest(3,:).*parest(2,:).^3); -.5./sqrt(parest(2,:).*parest(3,:).^3)];
    cn = [zeros(size(parest(1,:))); zeros(size(parest(1,:))); ones(size(parest(1,:)))];
    cc = sparseblock([cm(:), cs(:), cn(:)],3,'transpose');    
    S = - cc'*(D2L(parest))^-1*cc;
    err = full(reshape(sqrt(diag(S)),3,ncol));
    
end
%%% Change to usual parameterization
parest(2,:) = (parest(2,:).*parest(3,:)).^-.5;

if nargout > 2
    lik = L(parest);
end