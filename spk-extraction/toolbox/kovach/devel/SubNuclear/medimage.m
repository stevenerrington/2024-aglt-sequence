

classdef medimage <handle

    %% Medical image class 
    %
    % A class to load medical images and store related data.
    
% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/medimage.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

    properties
        
        Label ='';  %% Image label
        header = [];  %% Header data
        vox2mm = [];  %% Vox to mm affine trans. matrix
        vox2tal = [];  %% Vox to talairach or mni transfm. matrix.
        path='';  %%Path to file
        file='';  %% file name
        subjectID='';  %% subject/patient ID
        userdata = [];  %% Other data from user
        notes = '';
    end
    
    properties (Dependent=true)
     Data
     load2mem %% Load image data to memory if true, otherwise read from disk.
  
    end
    
    properties (Hidden=true)
        
        imgDat = [];
        inmem = true;
    end
    
    methods
        
        function me = medimage(varargin)
            
            global MEDIM_LASTPATH
            
            if ~isempty( MEDIM_LASTPATH)
                me.path = MEDIM_LASTPATH;
            end
            
            if nargin == 1 && isempty(varargin{1}) %%% return a volumeview without loading data if the only input is an empty matrix.
                me = me([]);
                return 
            end
            if nargin > 0
                me.Label = varargin{1};
            end
            if nargin > 1
                me.file = varargin{2};
            end
            if nargin > 2
                me.path = varargin{3};
            end
            if nargin <2
                me.asgnData;
            else
                me.asgnData(fullfile(me.path,me.file));
            end
                
            MEDIM_LASTPATH = me.path;
            
        end
        
        
        %%%%%%%%
        function asgnData(me,a)
            
            
            if nargin < 2 || isempty(a)
                
                 [fn,pth] = uigetfile({'*.nii*','*.hdr';'nifti','img'}',sprintf('Find a %i image for %s',me.Label, me.subjectID),me.path);
     
            else
                [pth,fn,ext]= fileparts(a);
                fn = [fn,ext];
                if nargin > 2 
                    me.path = pth;
                end
            end                
            me.path = pth;
            me.file = fn;
            me.loadfile;
        end
            
        function set.Data(me,a)
            asgnData(me,a);
        end
        
        %%%%%%%%%%% 
        function x = get.Data(me)
            
            if me.inmem
                x=me.imgDat;
            else
                [x,me.header] = loadfile(me);
            end
        end
        %%%%%%%
        function set.load2mem(me,a)
            if ~islogical(a) && ~ismember(a,[0 1])
                error('Input must be true (1) or false (0)')
            end
            me.inmem = a;
            if a
                [me.imgDat,me.header] = me.loadfile;
            else
                me.imgDat=[];
            end
        end
        
        %%%%%%%
        function a=get.load2mem(me)
            a= me.inmem;
        end
        %%%%%%%%%%% load image
        function [img,hdr] = loadfile(me)
              [p,f,ext ]= fileparts(me.file);
              fn = fullfile(me.path,me.file);
               switch lower(ext)
                    case {'.hdr'}
                        [img,hdr] = loadimg(fn);

                    case {'.nii'}
                          [img,hdr] = readnifti(fn);    
                   case {'.gz'}
                       gunzip(fn)
                       me.file =fullfile(p,f);
                       me.loadfile();
                       return
                    otherwise
                    error('Unrecognized file extension %s',ext)
               end
               perm = [1 2 3];
               e = eye(4); e = e([perm,4],:);
              img = permute(img,perm);
               me.header = hdr;
               me.vox2mm = e(1:3,1:3)*hdr.vox2unit(perm,:)*e'^-1;
               if me.load2mem
                   me.imgDat = img;
               end
               
        end   
        
            
    end
end
