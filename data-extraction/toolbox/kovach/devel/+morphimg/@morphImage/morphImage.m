classdef  morphImage < handle

% 
% IMMORPH morphs images by applying a piecewise linear transformation over a delaunay tesselation.
%   The image interpolation is carried out by GRIDDATA
%
% How to use:
% 
% 1. load two images (X and Y) into Matlab  using imread, eg. X = imread('imageName.ext','ext').
% 2. Convert images to double and scale to [0 1], eg  X = double(X)/255.
% 3. Add fiducial points with the immorph command, M = immorph(X,Y)
%         
%         A. Right clicking will add a new active point to each image. The active point is indicated with a cross.
%         B. Left clicking will move an active point to a new position.
%         C. Right clicking again will de-activate the current active point.
%         D. Left clicking over an inactive point will activate it if no other point is already active.
% 
% 4. Press the space bar when finished adding fiducial points.
% 5. IMMOPRH returns a data structure, M, with the following fields
%         
%             X: matrix with the original image X
%             Y: matrix with the original image Y
%             VX: fiducial point locations on X
%             VY: fiducial point locations on Y
%             Tess: matrix containing indices of the ponts in VX and VY that form the vertices of each tesselation.
%                 Tesselation is done with the DELAUNAY function.
%             mapto: mapping of each pixel in X to a pixel in Y
%             makepic: function handle to the function that generates the morph.
%             tessmap: function handle to a function that generates a map of tesselations.
% 
% 6. IMMORPH takes the following optional keywords
% 
%     'align': Aligns Y with X by rotating and scaling Y to minimize squared difference between the fiducial points.
%     'vertices': Followed by VX and VY. Adds vertex (fiducial) points from the matrices VX and VY:  
%         M = immorph(X,Y,'vertices',VX,VY)  adds points in VX to X and corresponding points VY to Y.These can be edited subsequently.
%     'noget': Runs IMMORPH without acquiring additional fiducial points. Uses only points supplied with the
%         'vertices' option.
%     'fix': Followed by 1 or 2. Disallows any changes to points in X ( if 1 ) or Y ( if 2 ). No additional points can be
%         added to either image if either of them is fixed.
%         
%     
% Generating Morphs:
% 
%    1. Morphs are generated with the M.makepic function handle, which is called using FEVAL:
%      The structure M, itself, has to be passed as an argument to M.makepic. 
% 
%         morphPic = feval(M.makepic,M,morphlevel).
%      
%       where morphlevel is the fractional change between X and Y (morphlevel = .5 gives a 50%
%           morph).
%
%     2. Morphing level can be specified independently for each pixel by using the keyword 'mask' followed by 
%        a matrix the same size as the image: 
%               
%        morphPic = feval(M.makepic,M,morphlevel,'mask',maskMatrix)
% 
%     3. Shading and warping can be controlled independently using the keywords 'shadingscale' and 'warpingscale' respectively.
%         In both cases 0 preserves X and 1 preserves Y. To warp X to Y use:
%         
%             morphPic = feval(M.makepic,M,0,'warpingscale',1)  or 
%             morphPic = feval(M.makepic,M,1,'shadingscale',0)

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/imagetools/immorph.m $
% $Revision: 89 $
% $Date: 2012-03-28 00:47:56 -0500 (Wed, 28 Mar 2012) $
% $Author: ckovach $
% ------------------------------------------------

% Written by Christopher Kovach 2006,2012
  


    properties ( SetAccess = private )        
%         VX = [];
%         VY = [];
%         Y = [];
%         X = [];
        images;
        Tri = [];
%         mapto = [];
        axes=[];
        trh = [];
        tesswgt = [];
        detA = [];
    end
    
    methods 
        function morphim = morphImage(varargin)
           
            createMorphMap(morphim,varargin{:})
            
        end
    end
            
    methods ( Access = private )
        
        [VX,VY] = getvertices(morphim,varargin)
        
        [Yout,VYout] = alignImages(morphim);
        
        trplot(bh,ev,hd)
    end
    
end
    


