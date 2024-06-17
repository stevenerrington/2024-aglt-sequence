function C = cscale(obj,varargin)

%C = indscale(obj) returns a matrix of indices of the colormap to which each corresponding element of the 
%unscaled data represented in object obj is mapped. If the keyword 'RGB' is included, then 
%C is the nXmX3 TrueColor rgb representation of axis ax. Obj must be an axis child and must have a 
%valid CData field.

%C. Kovach 2002


i = 1;
convert = 0;

varg = varargin;

while i <= length(varg)
  switch lower(varg{i})
    case 'rgb'
        convert = 1;
        i = i+1;
    otherwise
        error([varg{i},' is not a valid keyword'])
  end
end

ax = get(obj,'parent');
fig = get(ax,'parent');

X = get(obj,'CData');
clim = get(ax, 'clim');
cmap = get(fig, 'colormap');
csz = size(cmap,1);

X(find(X<clim(1))) = clim(1);
X(find(X>clim(2))) = clim(2);

C = X - clim(1);
C  =  round(C*(csz-1)/diff(clim)) + 1;

if convert
    C = reshape(cmap(C,:),size(C,1),size(C,2),3);
end