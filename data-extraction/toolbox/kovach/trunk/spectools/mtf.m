
function [T,C,shapemat,jkT,jkC] = mtf(X,varargin)

% function [T,C] = mtf(X,S1,S2,...SK)
% Computes transfer functions in the case of multiple inputs. 
% X is a TxN output, where each T is a window across which will be
% averaged. N must be at least equal to the number of inputs, and each 
% S must be a matrix of the same size as X representing an input.
% C is the estimate of the cross spectrum. 
%
% See also MTCSPECT

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/spectools/mtf.m $
% $Revision: 247 $
% $Date: 2013-07-23 17:20:13 -0500 (Tue, 23 Jul 2013) $
% $Author: ckovach $
% ------------------------------------------------

%C. Kovach 2008

win = 'hanning';

i = 1;
while i <= length(varargin) && isnumeric(varargin{i})
    inputs{i} = varargin{i};
    i = i+1;
end

noisespect = 1e-9;
winN = size(X,1);
noverlap = 0;
nw = 2;
multitaper = true;  % Currently this only works with multitaper
jacknife = false;
shapemat = [];
jkT = [];
jkC = [];
princomp = 0;
while i <= length(varargin)
    switch lower(varargin{i})
       case 'window'
            win = varargin{i+1};
            i = i+1;
            multitaper = false;
       case 'winn'
            winN = varargin{i+1};
            i = i+1;
       case 'noverlap'
            noverlap = varargin{i+1};
            i = i+1;
       case 'nw'        %Specifies bandwidth parameter for  DSPSS windows 
            nw = varargin{i+1};
            multitaper = 1;
            i = i+1;
       case 'noisespect'        %Eventually do this using multitaper DSPSS windows 
            noisespect = varargin{i+1};
            i = i+1;
       case 'multitaper'        
            multitaper = varargin{i+1};
            i = i+1;
       case 'princomp'        
            princomp = varargin{i+1};
            i = i+1;

     case 'jacknife'        
            jacknife = varargin{i+1};
            i = i+1;
%        case 'jacknife'        
%             jacknife = true;
%            
       otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end


if isempty(win) || strcmpi(win,'none')
    win = 'rect';
end

% 
% if princomp~=0
%     
%     if princomp>1
%    
% %         [U,~] = svds(fft(cat(1,X,zeros(size(X) - [1 0]))),princomp);
%         [U,~] = svds(fft(cat(1,X,X(end:-1:2,:))),princomp);
%     elseif princomp>0
%         [U,v] = svds(cat(1,X,zeros(size(X) - [1 0])),min(size(X)));
%         [U,v] = svds(cat(1,X,zeros(size(X) - [1 0])),min(size(X)));
%         U = U(:,cumsum(diag(v)) <= sum(diag(v))*princomp);
%     end
% %     U = fft(cat(1,U,zeros(size(U) - [1 0])));
% %     U = fft(U);
% else
    U = 1;
% end

ninpt = length(inputs); %number of inputs


if min(size(X)) ~= 1
    
    nwins = size(X,2); %number of windows
    nwinpts = size(X,1);
    reshapex = false;
else
    nwinpts = winN;

    winstart = 1:nwinpts-noverlap:length(X)-nwinpts+1;    
    nwins = length(winstart);
    winrg = (0:nwinpts-1)';
    
    shapemat = repmat(winstart,nwinpts,1) + repmat(winrg,1,nwins);
    X = X(shapemat);
    
    reshapex = true;

end
    
    
IN = zeros(size(X).*[1 ninpt]);
for i = 1:ninpt    
%     S =  diag(sparse(inputs{i}(:)));
%     S =  
    if reshapex
        
        IN(:,(1:ninpt:end) + (i-1)) = inputs{i}(shapemat);
    else
        IN(:,(1:ninpt:end) + (i-1)) = inputs{i};
    end        
end

%Window, pad and apply DFT
if multitaper 
    [W, weights] = dpss(winN,nw,2*nw-1);
elseif ischar(win)
    try
%       W = repmat(window(nwinpts,win),1,nwins);
      W = window(win,nwinpts);
    catch
%       W = repmat(window(win,nwinpts),1,nwins);
      W = window(win,nwinpts);
    end   
    weights = 1;
else
%       W = repmat(win,1,nwins);
      W = win;
      weights = ones(1,size(W,2));
end

XXs = spalloc(((size(X,1)*2-1)*ninpt)^2,1,0);
XYs =  spalloc(((size(X,1)*2-1))^2*ninpt,1,0);
% if jacknife
% jkXX = repmat({sparse(0)},1,nwins);
% jkXY = repmat({sparse(0)},1,nwins);
% %     jkXX = sparse([nxwftpts nxwftpts);
% %     jkXY= sparse([nxwftpts nxwftpts);
% end
%     
for i = 1:size(W,2)
    ww = repmat(W(:,i),1,nwins);
    
    INw = cat(1,IN.*repmat(ww,1,ninpt),zeros(size(IN) - [1 0]));
%     INw = IN.*repmat(ww,1,ninpt);
    
    INwft = U'*fft(INw);
    
    Xwft = U'*fft(cat(1,X.*ww,zeros(size(X) - [1 0])));
%     Xwft = fft(X.*ww);
    
    nxwftpts = size(Xwft,1);

    %Rearrange INft into a block diagonal matrix
    K1 = kron(ones(1,nwins),speye(nxwftpts*ninpt));
    K2 = kron(speye(nwins),kron(ones(ninpt,1),speye(nxwftpts)));
    K3 = kron(ones(1,nwins),speye(nxwftpts));
    
    

    spINft = K1*diag(sparse(INwft(:)))*K2;
    spXft =  K3*diag(sparse(Xwft(:)));
    
    spINfts{i} = spINft;
    spXfts{i} = spXft;
    %Compute the power and cross spectral power between inputs
% %     XX = spINft*spINft' + weights(i)*XX;
%     XX = spINft*spINft'*weights(i) + XX;
    m2lin = @(x)x(:);
    XXs(:,i) = m2lin(spINft*spINft'*weights(i));    
    %Compute cross spectrum between outputs and inputs
% %     XY = spXft*spINft' + weights(i)*XY;
%     XY = spXft*spINft'*weights(i) + XY;
    XYs(:,i) = m2lin(spXft*spINft'*weights(i));
% 
%     if jacknife
%        
%         spI = speye(nxwftpts*ninpt);
%             fprintf('\nComputing jacknife...\n')
%         for j = 1:nwins
% %             fprintf('%4i',j)
%             trind = zeros(1,nwins);
%             trind(j) = 1;
%             jkK1 = K1 - kron(trind,spI) ;
%             jkK2 = K2 - kron(diag(trind),kron(ones(ninpt,1),speye(nxwftpts)));
%             jkK3 = K3 - kron(trind,speye(nxwftpts));
%             
%             jkspINft = jkK1*diag(sparse(INwft(:)))*jkK2;
%             jkspXft =  jkK3*diag(sparse(Xwft(:)));
%             
%             jkXX{j} = jkspINft*jkspINft'*weights(i) + weights(i)*jkXX{j}; 
%             jkXY{j} = jkspXft*jkspINft'*weights(i) + jkXY{j};
%             
% %             trind(i) = 0;
% %             fprintf('\b\b\b\b')
%         end
%              fprintf('Done.\n')
%     end
end
XX = reshape(sum(XXs,2),size(K1,1)*[1 1]);
XY = reshape(sum(XYs,2),[size(K3,1) size(K1,1)]);

if length(noisespect) == 1 && noisespect ~=0
    noisespect = noisespect*speye(size(XX));
elseif  size(noisespect,2)== 1 
    noisespect = diag(sparse([noisespect;noisespect(1:end-1)]));
end

XX = XX + sparse(noisespect); %add noise related regularization 



%And Voila! The LSE estimate for the transfer function
%T = XX\XY'; Oops. Inverse is efficient for Cholesky factor

if princomp~=0  %Principle component dimensionality reduction
    
    if princomp>1
   
%         [U,~] = svds(fft(cat(1,X,zeros(size(X) - [1 0]))),princomp);
        [U,~] = svds(XX,princomp);
    elseif princomp>0
        [U,v] = svds(XX,min(size(X)));
        U = U(:,cumsum(diag(v)) <= sum(diag(v))*princomp);
    end
%     U = fft(cat(1,U,zeros(size(U) - [1 0])));
%     U = fft(U);
else
    U = 1;
end
forcepd= @(x)(x+x')/2; %Force matrix to be semidefinite
xx = forcepd(U'*XX*U);
R = chol(xx);
xy = XY*U;
B = (xy/R/R');
Tdg = B*U';

T = [];
for i = 1:ninpt   
    T(:,i) = diag(Tdg(:,(i-1)*nxwftpts+(1:nxwftpts)));
end

if jacknife
   

    spI = speye(nxwftpts*ninpt);
        fprintf('\nComputing jacknife...\n')
    jkT = zeros(size(T,1),nwins);
    for j = 1:nwins
        jkXX = 0;
        jkXY = 0;
 %             fprintf('%4i',j)
       for i = 1:size(W,2)
            trind = zeros(1,nwins);
            trind(j) = 1;
            jkK1 = K1 - kron(trind,spI) ;
            jkK2 = K2 - kron(diag(trind),kron(ones(ninpt,1),speye(nxwftpts)));
            jkK3 = K3 - kron(trind,speye(nxwftpts));

            jkspINft = jkK1*diag(sparse(INwft(:)))*jkK2;
            jkspXft =  jkK3*diag(sparse(Xwft(:)));
            
            Q = U'*jkspINft;
            jkXX = Q*Q'*weights(i) +jkXX; 
            jkXY = jkspXft*Q'*weights(i) + jkXY;
       end
%             trind(i) = 0;
%             fprintf('\b\b\b\b')
 
        R = chol(forcepd(U'*(jkXX+ sparse(noisespect))*U));
        Tdg = ((jkXY{jj}*U)/R/R')*U';
        for k = 1:ninpt   
             jkT(:,j,k) =   U*diag(Tdg(:,(k-1)*nxwftpts+(1:nxwftpts)));
        
            if nargout > 4
                for kk = 1:ninpt       
                    jkC(:,kk,k) = U*full(diag(jkXY{j}((kk-1)*nxwftpts+(1:nxwftpts))));
                end
            end

        end
    end
    
     fprintf('Done.\n')
end



if nargout > 1
    for i = 1:ninpt       
        C(:,i) = full(diag(XY(:,(i-1)*nxwftpts+(1:nxwftpts))));
    end
end


