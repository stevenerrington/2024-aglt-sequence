%%%%%%%%%%
function [pk,pksign,zeroc] = getpeak2(x)

% Extract all positive and negative peaks and zero crossings 

tol = eps;

dx = diff([x(end);x;x(1)]);  %Again, compute difference on circular domain

dx = dx.*(abs(dx)>tol); %Apply tolerance threshold

sdx = sign(dx) ;

dsdx = [sdx(end);sdx] - [sdx;sdx(1)];

zeroc = sign([x;x(1)])~=sign([x(end);x]);
% dsdx = diff(sdx);

pk = abs(dsdx)>0;
% pk([1 end]) = true;
pksign = sign(dsdx.*pk);

%%% When one or more adjacent values are equal, assign a peak  to
%%% the first inflection for concave or convex regions, and to both
%%% inflections for saddle points.
if any(dx==0) && any(dx~=0) 
    fpks = pksign(pk);
    pk(pk) = [fpks(end);fpks(1:end-1)] ~= fpks;
    pksign = pk.*pksign;
end
    
    

pk = pk(2:end-1);
pksign = pksign(2:end-1);
zeroc = zeroc(2:end-1);

