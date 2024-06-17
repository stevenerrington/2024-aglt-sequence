function clout = splitmogfit(X,hos,n,npc,plh)

if nargin <3
    n = 2;
end
if nargin < 4
    npc = ceil(length(hos)*1.5);
end
if nargin < 5
    plh = true;
end
min_clust_size = 10;

[A,~,~,~,~,plh]  = hos.get_block(X,25,plh);
if all(plh(1).YData==0)
    plh(1).YData = 1:size(X,2);
end
    
% XF = hos.xfilt(A(:,:,1));
win = exp(-fftshift(hos(1).sampt/hos(1).sampling_rate/(1e-3)).^2/2);
XF = fftshift(hos.xfilt(X),1).*win';
XF = XF./sqrt(nanmean(nanmean(XF.^2,1),2)+eps);
XF = reshape(permute(XF,[1 3 2]),size(XF,1)*size(XF,3),size(XF,2));
[~,~,v] = svds(XF,npc);

%  [cl,gm,gms] =   mogclust(v,1:n,'Replicates',10,'Options',struct('MaxIter',500));
 cl = isosplit5(v',struct('isocut_threshold',.4));
 
% c2c =   chi2_criterion(v,1:n);
%c2c =   chi2_criterion(v,1:n,@(x,k)kmeans(x,k,'Replicates',20,'Options',struct('MaxIter',500)));
% cl = c2c.OptimalY;
if max(cl)==1
    clout = cl;
    return
end
clout = zeros(size(X,2),1);

srt = sort(cl);
 yd = plh(1).YData;
for k = 1:max(cl)
    try
        if sum(cl==k)<min_clust_size
            error('Minimum cluster size reached')
        end

         plh(1).YData = yd(srt==k);
        cl2 = splitmogfit(X(:,cl==k,:),hos,n,npc,plh);
    catch err
        warning(err.message)
        cl2 = 1;
    end
%     clout(cl==k) = k + (cl2-1)/n;
    clout(cl==k) = k + (cl2-1)/max(cl);
    
end
    