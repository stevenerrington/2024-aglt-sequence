function [WIN,Y] = MakeMask(X, varargin)

%[W,Y] = MakeMask(X);
%
% A semiautomated script for masking image parts.
%
% X is the image matrix
%
% W is the mask with ones in the foreground and zeros in the background
%
% Y foreground of X with a gray background.
%
% This function detects edges as the inverse of the Laplacian weighted by the
% gradient slope, ie  sqrt((dX/dx).^2 + (dX/dy).^2) ) / del2( X ) > threshold.
%
% The user repairs gaps in the edges, then identifies regions of the thresholded map 
% which are part of the object of interest. The mask is computed from the
% convex hull of the edges of the area of interest (only convex masks can be computed
% by this script). 
%

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/imagetools/MakeMask.m $
% $Revision: 39 $
% $Date: 2011-06-04 23:21:09 -0500 (Sat, 04 Jun 2011) $
% $Author: ckovach $
% ------------------------------------------------

%C Kovach 2006

% Main paramaters affecting the performance of the algorithm
threshlo = 0;  %Thresholds for the high smoothed image
threshhi = 2e1;  % Threshold for the low smoothed image
klo = 6; %Smoothing kernel width in pixels for high smoothing
khi = 6; %For lower smoothing

i = 1;
while i < length(varargin) 
    
    switch lower(varargin{i})
        
        case 'threshlo'
            threshlo = varargin{i+1};
            i = i+1;
        case 'threshhi'
            threshhi = varargin{i+1};
            i = i+1;
        case 'smoothinglo'
            klo = varargin{i+1};
            i = i+1;
        case 'smoothinghi'
            khi = varargin{i+1};
            i = i+1;
        otherwise
            error(sprintf('%s is not a valid keyword',varargin{i}));
    end
    i = i+1;
end
        
    

X = double(mean(X,3)); %in case the image is truecolor 

%Gaussian smoothing kernel
glo = normpdf(-4*klo:1:4*klo,0,klo);
Glo = glo'*glo;

ghi = normpdf(-4*khi:1:4*khi,0,khi);
Ghi = ghi'*ghi;



Xsm = conv2(X,Glo,'same'); %smoothed image
Xhi = conv2(X,Ghi,'same'); %Higher frequency

[gx,gy] = gradient(Xsm);
slope = sqrt(gx.^2+gy.^2);

[gxhi,gyhi] = gradient(Xhi);
slopehi = sqrt(gxhi.^2+gyhi.^2);


%Laplacian
D2sm = del2(Xsm);
D2hi = del2(Xhi);

%Convex side of zero crossing taken as gradient steepness divided by del2
ZC = (slope./(D2sm) > threshlo).*(slopehi./D2hi > threshhi);

ZC = draw(ZC); %User fills gaps

[Edges,F] = fill(ZC);

%Get Convex hull from edges and inside part

[I,J] = ind2sub(size(Edges),find(Edges));
k = convhull(I,J);

vi = I(k);
vj = J(k);

I = [1:size(Edges,1)]'*ones(1,size(Edges,2));
J = ones(size(Edges,1),1)*[1:size(Edges,2)];

WIN = inpolygon(I,J,vi,vj);

Y = X;
Y(~WIN) = 133;

imagesc(Y),axis image, colormap gray
figure(gcf)


%%%%%%%%%%%%%%%%%%%%%%
function [Edges,F] = fill( ZC )

%A function for filling and finding the edges of regions of continuous zeros in a matrix
figure
F = zeros(size(ZC));
ZC(:,1) = 1; ZC(:,end) = 1;
ZC(1,:) = 1; ZC(end,:) = 1;

while 1

    imagesc(ZC+2.*F), axis image, axis ij, caxis([0 2])
    title('Choose a point in the foreground. Right-click to finish.')
    pt = ginput(1);
    if strcmp(get(gcf,'selectiontype'),'alt')
        break
    end
    open = 1;

    F(round(pt(2)),round(pt(1))) = 1;

    Edges = zeros(size(F));

    while open

        %% I direction
        for perm = 1:2

            permat = [2 1];

            F = permute(F,permat);
            ZC = permute(ZC,permat);
            Edges = permute(Edges,permat);

            Q = zeros(prod(size(ZC)),1);

            filldi = find(diff(F(:),1));

            zcdi = find(diff(ZC(:),1));

            %samez = ismember(zcdi,filldi);
            samef = ismember(filldi,zcdi);
            Edges(filldi(samef) + 1) = 1; 

            %zcdi(samez) = [];
            filldi(samef) = [];

            if isempty(filldi)
                open = 0;
            else
                open = 1;
            end


            spots = cat(1,cat(2,filldi,ones(size(filldi))),cat(2,zcdi,zeros(size(zcdi))));

            [q, sortind] = sort(spots(:,1));

            sorted = spots(sortind,:);

            holes = diff(sorted(:,2));

            Q(sorted(find(holes == 1),1)+1) = 1;
            Q(sorted(find(holes == -1)+1,1)+1) = -1;
            tofill = logical(cumsum(Q));

            F(tofill) = 1;


        end
    end

end

close(gcf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = draw(ZC)


global Z

Z = ZC;

figure, imagesc(Z), axis image, colormap gray
title('Repair gaps. Right click to finish.')
axis tight


set(gcf,'windowbuttonupfcn',@bup);
set(gcf,'windowbuttondownfcn',@bdown);


uiwait(gcf)

%%%%%%%%%%%%
function x = setpoint(src,evnt,ZC)

% Fills points based on pointer position

pensize = 2; %1/2 width of the pen.

global Z

psel = [1 0]*get(gca,'currentpoint')*[0 1; 1 0; 0 0];
    %psel(:,2) = get(gca,'currentpoint')*[1 0; 0 1; 0 0];
psel = round(psel);
rg1 = psel - pensize;
rg2 = psel + pensize;

rg1 = rg1.*(rg1 > 0) + psel.*(rg1 <= 0);
rg2 = rg2.*(rg2 <= size(Z)) + psel.*(rg2 >(size(Z)));

if psel > 0 & psel <= size(Z);
    Z(rg1(1):rg2(1),rg1(2):rg2(2)) = 1;
end

imagesc(Z), axis image, colormap gray
drawnow

%%%%%%%
function donothing(src,evnt)

return

%%%%%%%
function bdown(src,evnt)

%Callback for buttonpress over the drawing figure
global Z

if strcmp(get(gcf,'selectiontype'),'alt')
  delete(src)
else
 set(gcf,'windowbuttonmotionfcn',{@setpoint,Z})
 setpoint(src,evnt,Z)
end


%%%
function bup(src,evnt)

set(gcf,'windowbuttonmotionfcn',@donothing)




