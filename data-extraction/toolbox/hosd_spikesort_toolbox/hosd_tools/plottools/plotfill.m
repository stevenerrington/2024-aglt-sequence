function varargout = plotfill(x,y,c)

%creates a filled plot

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 2 
    y=x;
    x = [];
end

if min(size(y))==1
    y = y(:);
end
if isempty(x)
    x = (1:size(y,1))';
end
if min(size(x))==1
    x = x(:);
end

if nargin <3
    c = [0 0 1; 0 1 0; 1 0 0; 1 1 0; 1 0 1; 0 1 1];
end

if min(size(x))==1 && min(size(x))~=min(size(y));
    x = repmat(x(:),1,size(y,2));
end

x = x([1, 1:end,end],:);
y = cat(1,zeros(1,size(y,2)),y,zeros(1,size(y,2)));
hold on
for i = 1:size(y,2)
    h(i) = fill(x(:,i),y(:,i),c(mod(i-1,size(c,1))+1,:));
end

if nargout > 0
    varargout{1} = h; 
end
