
function [ptch,hgt] = connectmap(pos,C,thresh,scaling,starthgt,scalewidth)

% function [ptch,hgt] = connectmap(pos,C,thresh,scaling,starthgt,scalewidth)
% 
%  Plots connections between points at locatons in 'pos', contained in matrix 
% 'C' as bars scaled in width and color to the values in C. 
%
% input arguments:
% 
% pos    -  N x 2 matrix of locations
% C      -  N x N matrix of connection weights
% thresh -  Either (1) a scalar threshold such that only connections for which
%           abs(C(i,j)) > thresh are plotted, or (2) an N x N binary matrix 
%           such that only connections for which thresh(i,j) == 1 are
%           plotted (default 0).
% scaling- scaling factor between bar width in pixels and the values in C (default .1).
% starthgt- minimum height in the z plane at which to plot bars (higher bars lie on
%            top of lower bars (default 0).
% scalewidth - bar width is not scaled to C(i,j) if false, but is constant (default true).
% 
% 
% output arguments:
% 
% ptch - handles for patch objects for each bar
% hgt  - maximum height of bars in the z-plane.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


% C. Kovach 2009


if nargin < 6
    scalewidth = true;
end

if nargin < 5
    starthgt = 0;
end
if nargin < 3
    thresh = 0;
end

if nargin < 4
    scaling = .1;
end
ndim = size(pos,2);
ptch = [];
srtc = max(abs(C),abs(C'));  %sort according to absolute value
srtc = srtc - diag(diag(srtc));
[srcoh,csrt] = sort(srtc(:));

[J,I] = meshgrid(1:size(C,1),1:size(C,2));

hgt = starthgt;
for i = find(srcoh)'
        p1 = I(csrt(i));
        p2 = J(csrt(i));
        
        cs = [C(p1,p2),C(p2,p1)];
        if length(thresh) == 1 && max(abs(cs)) < thresh || p1==p2  % fixed threshold
            continue
        elseif length(thresh)>1 && ~thresh(p1,p2)  % significance defined for each pair
            continue
        end
        
        dp = diff(pos([p1 p2],:));
        dpn = dp./sqrt(sum(dp.^2));
        orthmat = eye(ndim) - dpn'*dpn; 
        [u,e] = eigs(orthmat);
        U = u(:,diag(e)>1e-10);
        
        if scalewidth
            CC = C;
        else 
            CC = ones(size(C));
        end
        L1 = [U*abs(CC(p1,p2)).*scaling + repmat(pos(p1,:)',size(U,2)),-U*abs(CC(p1,p2)).*scaling + repmat(pos(p1,:)',size(U,2))];
        L2 = [U*abs(CC(p2,p1)).*scaling + repmat(pos(p2,:)',size(U,2)),-U*abs(CC(p2,p1)).*scaling + repmat(pos(p2,:)',size(U,2))];
        
        for k = 1:2*size(U,2)
            ll = [L1(:,[k mod(k,size(L1,2))+1]),L2(:,[mod(k,size(L2,2))+1 k])];
            ll(4,:) = 0;
            cols = [cs(1)*ones(1,2*size(U,2)),cs(2)*ones(1,2*size(U,2))];
            
            ptch(end+1) = patch( ll(1,:),ll(2,:),ll(3,:) + hgt,cols);
%             shading interp
            if size(U,2) == 1
                hgt = hgt+.001;%10*eps;
            end
        end
    
end

        
        
        
        
        
        