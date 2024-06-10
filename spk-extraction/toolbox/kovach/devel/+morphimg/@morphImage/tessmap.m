function [TM,Ind] = tessmap(morphim,pcnt)

% A function that returns a map of pixels showing the index of the
% tesselation to which they belong.

XX = cat(4,morphim.images.X);
VX = cat(3,morphim.images.vert);
szeX = size(XX);
% szeY = size(morphim.Y);

Tess = morphim.Tri;

% VX = morphim.VX;
% VY = morphim.VY;

Ind = logical(spalloc(prod(szeX(1:2)),length(Tess),prod(szeX(1:2))));

MDX = VX(:,:,1)*(1-sum(pcnt));

nim = length(morphim.images);
for i = 2:nim
    MDX = MDX + pcnt(i-1)*VX(:,:,i);
end

[pixmapXI,pixmapXJ] = ind2sub(szeX(1:2),1:prod(szeX(1:2)));
% pixmapX = [pixmapXI;pixmapXJ];

TM = zeros(szeX(1:2));



for i = 1:size(Tess,1)
    


    insideptsX = inpolygon(pixmapXJ,pixmapXI,MDX(Tess(i,:),1),MDX(Tess(i,:),2));
    TM(insideptsX) = i;    
    Ind(:,i) = insideptsX;
  
    
end

