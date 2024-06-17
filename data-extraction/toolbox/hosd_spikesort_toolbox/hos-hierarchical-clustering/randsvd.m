
function [u,l,v] = randsvd(X,K)


sz = size(X);


reseed
rp = randperm(size(X,1));

u=zeros(K,size(X,2));
l=zeros(K);

n = 0;
xx = 0;
for k = 1:K:length(rp)
    
    x = X(rp(k:min(end,k+K-1)),:);
    
   xx=x*x';

    [vx,lx] = svd(xx);
    
    
   
    ux = diag(sqrt(diag(lx)).^-1)*vx'*x;
    
    if k >1
       ux = diag(sign(diag(u*ux')))*ux;
    end
    u(1:size(x,1),:) = (u(1:size(x,1),:)*n+ux(1:size(x,1),:))/(n+1);

    
     l(1:size(x,1),1:size(x,1)) = (l(1:size(x,1),1:size(x,1))*n+sqrt(lx(1:size(x,1),1:size(x,1))))/(n+1);

    n = n+1;
    k
end

% l = diag(sqrt(sum(u.^2,2)));

% u = u./sqrt(sum(u.^2,2));












































































































    
    
    
    
    
    