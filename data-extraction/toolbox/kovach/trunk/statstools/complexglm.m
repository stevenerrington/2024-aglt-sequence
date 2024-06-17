function [b,dev,pval,iXX,sigma,res,Yfit] = complexglm(Y,X,varargin)

% [b,dev,pval] = complexglm(Y,X,varargin)
% Fits gaussian GLM for complex valued data, returning coefficients,
% deviance and p-value compared to an intercept only model.
%
% If X is a 3-d array, then each column of X is regressed onto the
% corresponding column of Y, where the third dimension of X contains the
% independent measures.
%

i = 1;
intercept = true; % add intercept term if true
diagonly = false;
blockdiag = false;
reg = eps;
while i < length(varargin)
   switch lower(varargin{i})
       
       case 'intercept'
           intercept = varargin{i+1};
           i = i+1;
       case 'diagonly'  % Univariate within bands
           diagonly = varargin{i+1};
           i = i+1;
       case {'multiple','blockdiag'}
           blockdiag = varargin{i+1};
           i = i+1;
       otherwise
           error('Unknown keyword')
   end
    i = i+1;
end

N = size(X,1);
if size(X,3)>1 || blockdiag
    blockdiag = true;
    if intercept
        X(:,:,end+1) = 1;
        X0 = X(:,:,end);
        X0 = reshape(permute(X0,[1 3 2]),N,size(X,2)*size(X,3));
        X0 = sparseblock(X0,size(X,3));

    end
    nreg = size(X,3);
    X = reshape(permute(X,[1 3 2]),N,size(X,2)*size(X,3));
    X = sparseblock(X,nreg);
    Y = sparseblock(Y,1);
else
    if intercept
       X0=ones(N,1);
       X = [X,X0];
    end
    if ~diagonly
        nreg = size(X,2);
    else
        nreg =1;
    end
    
end

if intercept
    bintcpt = X0\Y;
    resintcpt = (Y-ones(N,1)*bintcpt);
    Gmod = sum(abs(resintcpt).^2);
    Cmod = sum(resintcpt.*resintcpt);
    Pmod = (abs(Gmod).^2 - abs(Cmod).^2)./abs(Gmod);
%     X = [X,ones(N,1)];
else
    resintcpt = Y;
    bintcpt =[];
    Gmod = sum(abs(resintcpt).^2);
    Cmod = sum(resintcpt.*resintcpt);
    Pmod = (abs(Gmod).^2 - abs(Cmod).^2)./abs(Gmod);
end

if ~diagonly
      if blockdiag

        XX = sparse(X'*X );
        XY = sparse(X'*Y);
        regmat = speye(size(X,2)*[1 1])*reg;


    %     iXX = inv(XX+regmat);
         chXX = chol(XX+regmat);
        iXX = chXX\(chXX'\speye(size(X,2)));

    %     b = chXX\(chXX'\(XY));
    %     b=X\Y;
        b = XX*(iXX*iXX)*XY;
        %%% Need to make sure that b is a sparse matrix with the right
        %%% structure for unsparsify to work. Unsparsify may return incorrect
        %%% values if called with the improper structure.
        nzindx = sparseblock(ones(nreg,size(Y,2)),1);  
         b(nzindx & b==0) = eps;
      else
          b = X\Y;
        if nargout > 3
            iXX = inv(X'*X);
        end
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
    dev = N*(log(Gmod)-log(Gres)+log(Pmod+eps)-log(Pres+eps));        
end
if blockdiag
    b = unsparsify(b);
    if intercept
        bintcpt = unsparsify(bintcpt);
    end
end
if nargout > 2
    pval = 1-chi2cdf(dev,(2-~any(imag(Y)))*(size(b,1)-size(bintcpt,1)));
end
if nargout > 4
    sigma = Gres./(N-1);
end
