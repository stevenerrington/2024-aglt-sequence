function [inlayer,onlayer] = insidebem(X,bem,outsidept)

% [inlayer,onlayer] = insidebem(X,bem,outsidept)
%Determines which boundaries contain points that are rows of X.
% It works by determining how many times the line between X and outsidept
% intersects each boundary. 

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/trunk/BEMtools/insidebem.m $
% $Revision: 186 $
% $Date: 2013-04-11 12:44:33 -0500 (Thu, 11 Apr 2013) $
% $Author: ckovach $
% ------------------------------------------------

%CKovach

if nargin < 3 
    outsidept = [100 100 100];
end

if isfield(bem,'layer')
    unq = unique(bem.layer);
else
    unq = 1:length(layer);
end

for i = 1: length(unq)
    
    if isfield(bem,'layer')
        if iscell(bem.pnt)
            pnt = [bem.pnt{:}];
        else
            pnt = bem.pnt;
        end
        tri = bem.tri(bem.layer == unq(i),:);
        [I,O] = intersects(X,repmat(outsidept,size(X,1),1),tri,pnt);
        
    else
        if iscell(bem(i).pnt)
            pnt = [bem(i).pnt{:}];
        else
            pnt = bem(i).pnt;
        end
        [I,O] = intersects(X,repmat(outsidept,size(X,1),1),bem(i).tri,pnt);
    end
    
    N = sum(I);
    
    inlayer(i,:) =  mod(N,2);
    onlayer(i,:) = sum(O)>0;
    
end

        
        
        