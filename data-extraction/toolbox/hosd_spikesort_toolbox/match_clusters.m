
function [rmapout,cmapout] = match_clusters(confmat,rfix,cfix)

% Match clusters by permuting the rows and columns of a confusion matrix
% until it is maximally diagonal.

if nargin < 2
    rfix = 1; %Leave this row unchanged
end
if nargin < 3
    cfix = 1; %Leave this column unchanged
end
sz = size(confmat);

rmap0=1:sz(1);
cmap0=1:sz(2);
rmap = setdiff(rmap0,rfix);
cmap = setdiff(cmap0,cfix);

r = exp(1i*(0:length(rmap)-1)/length(rmap)*2*pi-1i*pi);
c = exp(1i*(0:length(cmap)-1)'./length(cmap)*2*pi-1i*pi);

% figure
% im = imagesc(confmat);
 d = Inf;
 maxiter = 10;
 iter = 0;
 while d>0 && iter < maxiter
    
     newconfmat = confmat(rmap,cmap);
    m1 = r*newconfmat;    
    [srt,cmapi] = sort(angle(m1*conj(max(m1))));
    newconfmat = newconfmat(:,cmapi);
    
     m2 = newconfmat*c;    
     [srt,rmapi] = sort(angle(m2*conj(max(m2))));
     newconfmat = newconfmat(rmapi,:);
  
     d = norm(newconfmat-newconfmat(rmapi,cmapi));
    
     rmap = rmap(rmapi);
    cmap = cmap(cmapi);

 
%     im.CData = newconfmat;
%     drawnow
iter = iter+1;
 end

% %%% Shift to main diagonal
% sh = r*confmat(rmap,cmap)*conj(c)*exp(1i*pi/length(c));
% [~,cmapi2] = sort(angle(c.*sh));
% cmap = cmap(cmapi2);

[rmapi3,cmapi3] = sort_diag(confmat(rmap,cmap));
rmap = rmap(rmapi3);
cmap = cmap(cmapi3);


rmapout = rmap0;
rmapout(ismember(rmap0,rmap))=rmap;
cmapout = cmap0;
cmapout(ismember(cmap0,cmap))=cmap;


function [rmap0,cmap0] = sort_diag(confmat,rfix,cfix)

% Sort onto the diagonal

if nargin < 2
    rfix = []; %Leave this row unchanged
end
if nargin < 3
    cfix = []; %Leave this column unchanged
end
sz = size(confmat);

rmap0=1:sz(1);
cmap0=1:sz(2);
rmap = setdiff(rmap0,rfix);
cmap = setdiff(cmap0,cfix);


if isempty(rmap) && isempty(cmap)
    return
elseif isempty(rmap)
    [~,mxi] =max(confmat(1,cmap0));
    cmap0 = circshift(cmap0,1-mxi);
    return
elseif isempty(cmap)
    [~,mxi] =max(confmat(rmap0,1));
    rmap0 = circshift(rmap0,1-mxi);
    return
end

[mx,mxi] = max(confmat(rmap,cmap));
[~,mxi2] = max(mx);

rmap = circshift(rmap,1-mxi(mxi2));
cmap = circshift(cmap, 1-mxi2);

[rmap2,cmap2] = sort_diag(confmat(rmap,cmap),1,1);

rmap0(ismember(rmap0,rmap)) = rmap(rmap2);
cmap0(ismember(cmap0,cmap)) = cmap(cmap2);
