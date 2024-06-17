

function [m,k,stats] = vonmisesfit(X,unit)

%  [m,k] = vonmisesfit(X)
%
%   Maximum likelihood fitting of a von Mises (circular Gaussian)
%   distribution to the data in matrix X. Returns parameters m for mode and
%   k for dispersion.  X is assumed to be an angle in radians over [0,2pi].
%
%   This function uses an efficient vectorized method to fit data in each
%   column of X.
%    
%
%   [m,k,stats] = vonmisesfit(X)
%
%   Returns a structure with the p-value, likelihood ratio statistic, degrees of freedom.
%
%
%   [m,k,...] = vonmisesfit(X,data_type)
%
%   data_type = 'rad' Assumes X ranges over [0,2pi] (default)   
%   data_type = 'unit' Assumes X ranges over [0,1]   
%   data_type = 'phasor'  Assumes X is a complex phasor. The distribution 
%   is fitted to the phase angle (amplitude is ignored).
%  
%   If data_type is a scalar, then the data are multiplied by
%   2*pi/data_type.
%
%   NaN values in X are ignored.
%
%  See also RAYLFIT

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/statstools/vonmisesfit.m $
% $Revision: 56 $
% $Date: 2011-09-28 12:41:31 -0500 (Wed, 28 Sep 2011) $
% $Author: ckovach $
% ------------------------------------------------

% C. Kovach 2011


if nargin < 2  || isempty(unit)
    if any(abs(imag(X(:)))>0)
        unit = 'phasor';
    else
        unit = 'rad';
    end
end
    


switch unit
    case {'rad','radian','radians'}
        % Do nothing
    case {'unit'}
        X = 2*pi*X;
    case {'phasor'}  % convert a phasor to a phase angle
        PH = X./abs(X);
        X = atan2(real(PH),imag(PH));
    otherwise
        if ischar(unit)
            error('Unrecognized data type.')
        else
            X = X * 2*pi / unit;
        end
end

% X = X(~isnan(X));  %nans are ignored
tol = 1e-6;
Ns = sum(~isnan(X))';
ncol = size(X,2);

% trg = [cos(X(:)),sin(X(:))];
trg = zeros(size(X).*[1 2]);
trg(:,1:2:end) = cos(X);
trg(:,2:2:end) = sin(X);

trg(isnan(trg)) = 0;

% For the sake of numerical stability, we will estimate the canonical
% parameters for the exponential family,  a*cos(x) + b*cos(x), and
% compute the dispersion and mode from those.


if ncol == 1  % No need for sparse matrix with vector input
    
     kappas = @(x) sqrt(sum(x.^2)+eps);
     DM = trg;   % Design matrix with the canonical sufficient statistics for the exponential family
   
else
    if ~exist('sparseblock.m','file')
        error('For matrix inputs, this function requires the sparseblock function from the GazeReader Toolbox.\Inputs must be column vectors otherwise.')
    end
    ksum = kron(speye(ncol,ncol),ones(1,2));
    kappas = @(x) sqrt( ksum*x.^2 + eps );
    DM = sparseblock(trg,2);  
end

%%% Likelihood function
LLfun = @(pars) nansum(DM*pars ) - Ns'*( log(besseli(0,kappas(pars))) + log(2*pi) );  

%%% Derviatives of the bessel functions
DI = @(k) besseli(-1,k)./besseli(0,k);
D2I = @(k) ( besseli( -2 ,k ) + 1./(k+eps).*besseli(-1,k) )./besseli(0,k) -  (besseli(-1,k)./besseli(0,k)).^2; 

%%% Derivative of the likelihood function
LLdfun = @(pars) sum(DM)' - kron(Ns.*DI( kappas(pars) )./kappas( pars ),[1;1]).*pars;
       
%%% Second derivative of the likelihood function

if ncol == 1 % No need for sparse matrix with vector input

    LLd2fun = @(pars)   - Ns * ( D2I( kappas(pars) ).* pars*pars'/kappas(pars).^2 +...
                            DI( kappas(pars) ).* ( eye(2)./kappas(pars) - pars*pars'./kappas( pars ).^3 ) );

else
    
     spr = @(x)sparseblock(x',2)';

    Wfun = @(pars) sparse( Ns .* ( D2I( kappas(pars) )./kappas(pars).^2  - DI( kappas(pars) )./kappas( pars ).^3));
    LLd2fun = @(pars)   - spr(pars) * diag(Wfun(pars)) * spr(pars)' -...
                    diag(kron( DI( kappas(pars) ).*Ns./kappas(pars), sparse([1;1])  ) );
end                        
                            


pars = zeros(2*ncol,1);

dstep = Inf;

maxiter = 500;

iter = 0;

%%% Newton-Raphson gradient ascent.
while dstep > tol && iter < maxiter 
    

   del =  - LLd2fun( pars )^-1 * LLdfun( pars );
   pars = pars + del;
      
   iter = iter+1;
      
   dstep = sqrt(sum(del.^2));
   
   %    ll(iter)  = LLfun(pars);

end

if iter == maxiter
    warning('Von misses fitting reached maximum number of iterations (%i) and may mot be accurate',maxiter)
end

%%% This is the parameter, mu, giving the mode 
m = atan2(pars(2:2:end),pars(1:2:end));

%%% THis is kappa, giving the degree of phase dispesion (higher value ->
%%% lower dispersion)
k = kappas( pars );

if nargout > 2
    
    
    %%% Parameter estimates
    stats.efpars =  pars ;
    stats.I =   LLd2fun( pars ) ;
    stats.mkpars =  [m k] ;

    %%% Log Likelihood for each column of X
    
    LLindiv = @(pars) nansum( reshape( DM*pars,size(X,1),ncol )) - Ns'.*( log(besseli(0,kappas(pars))) + log(2*pi) )';  
    stats.LL =  LLindiv( pars ) ;
    stats.N =  Ns' ;

    %%% Likelihood ratio statistic
    stats.DDev = 2*(  LLindiv( pars )  -    LLindiv( zeros(size(pars)) ) );

    %%% Chi-squared degrees of freedom
    stats.df = 2;

    %%% P-value for null hypothesis (k = 0).
    stats.pval = 1 - chi2cdf(stats.DDev,2);

end
    
    