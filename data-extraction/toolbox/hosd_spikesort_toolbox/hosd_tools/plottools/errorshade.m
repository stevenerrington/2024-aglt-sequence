function [plh,h] = errorshade(x,y,stderr,varargin)


resamp2 = 100;
plotse = [0:1:2 2.5];
dist = @(x)exp(-x.^2/2);
basecolor = [1 1 1];
if length(y)<resamp2
    resamp2 = size(y,1);
end

plotwith='patch';

[a,b] = rat(size(y,1)/resamp2,.5);
se = resample(stderr,b,a);
xrs = x(round(1:a/b:end));

if strcmpi(plotwith,'patch')
     xpatch = [xrs;xrs(end:-1:1)];
     ypatch = [-ones(size(xrs));ones(size(xrs))];
end

plh = plot(x,y,varargin{:}); 
hold on
plotse = sort(plotse,'descend');
X = repmat(xrs,1,size(plotse,2));
for k = 1:length(plh)
    
    col = plh(k).Color;
    %ses = [plotse];
%     w = dist((plotse([1, 1:end-1])+plotse)/2);
      w = dist(plotse);
      
  
    yrs = resample(y,b,a);
    switch plotwith
        case 'pcolor'
            alpha = .5*ones(size(xrs,1),1)*w;
            Y = yrs*ones(1,length(plotse)) + se*plotse;
            C = repmat(permute(w'*col + (1-w)'*basecolor,[3 1 2]),size(xrs,1),1);

            h=pcolor(X,Y,Y);
            set(h,'FaceColor',col,'AlphaData',alpha*100,'FaceAlpha','interp','AlphaDataMapping','direct','EdgeColor','none');
        case 'patch'   
            alpha = .75*ones(2*size(xrs,1),1)*w;
            h = patch(xpatch(:,ones(1,length(plotse))),...
                        [yrs;yrs(end:-1:1)]*ones(1,length(plotse)) + ([se(:,k);se(end:-1:1,k)].*ypatch)*plotse,...
                        permute(col',[3 2 1]));
%                         permute(col'*w + basecolor'*(1-w),[3 2 1]));
            h.FaceVertexAlphaData = 100*alpha(:);
            h.FaceAlpha = 'flat';
            h.AlphaDataMapping = 'direct';
            h.EdgeColor = 'none';
    end
end
