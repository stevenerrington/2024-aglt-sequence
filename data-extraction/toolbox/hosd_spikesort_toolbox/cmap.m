
function cols = cmap(n)

%A colormap like hsv but spaced so the colors are more perceptually
%distinct.
cola = [0    0.2000    0.4000    0.6000    0.8000 1];
colb = [0    0.17    0.50    0.675    0.84 1];
colhsv = interp1(cola,colb,linspace(0,1-1/n,n+1));
cols = hsv2rgb([colhsv(1:n)',ones(n,2)]);

