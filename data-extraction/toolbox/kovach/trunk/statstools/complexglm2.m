function [b,dev,pval,iXX,sigma,res,Yfit] = complexglm(Y,X,varargin)

% [b,dev,pval] = complexglm(Y,X,varargin)
% Fits gaussian GLM for complex valued data, returning coefficients,
% deviance and p-value compared to an intercept only model.
%

i = 1;
intercept = true; % add intercept term if true
diagonly = false;
while i < length(varargin)
   switch lower(varargin{i})
       
       case 'intercept'
           intercept = varargin{i+1};
           i = i+1;
       case 'diagonly'  % Univariate within bands
           diagonly = varargin{i+1};
           i = i+1;
       otherwise
           error('Unknown keyword')
   end
    i = i+1;
end


if intercept
    bintcpt = ones(size(X,1),1)\Y;
    resintcpt = (Y-ones(size(X,1),1)*bintcpt);
    Gmod = sum(abs(resintcpt).^2);
    Cmod = sum(resintcpt.*resintcpt);
    Pmod = (abs(Gmod).^2 - abs(Cmod).^2)./abs(Gmod);
    X = [X,ones(size(X,1),1)];
else
    Gmod = 1;
    Pmod = 1;
    bintcpt = [];
end

if ~diagonly
    b = X\Y;
    if nargout > 3
    iXX = inv(X'*X);
    end
else
    iXX = diag(sum(abs(X).^2).^-1);
    b = diag(sum(Y.*conj(X))*iXX);
end

Yfit = X*b;
res = (Y-Yfit);
Gres = sum(abs(res).^2); % Covariance matrix
Cres = sum(res.*res); % "relation" matrix 
Pres = (abs(Gres).^2 - abs(Cres).^2)./abs(Gres);

%tsq = abs(b).^2./berr;

if nargout >1
%%% If interept is true then compute deviance as difference between full
%%% and intercept-only model, otherwise compute deviance as full model - 0 
    dev = size(X,1)*(log(Gmod)-log(Gres)+log(Pmod+eps)-log(Pres+eps));        
end
if nargout > 2
    pval = 1-chi2cdf(dev,(2-~any(imag(Y)))*(size(b,1)-size(bintcpt,1)));
end
if nargout > 4
    sigma = Gres./(size(X,1)-1);
end