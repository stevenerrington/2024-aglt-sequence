function pkmaps = find_contacts(template)
% function pkmaps = find_contacts(fname,template,sdthresh,sm)

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/ElectrodeIDtoolbox/find_contacts.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------


global pkmaps

[fname,path] = uigetfile({'*.nii','*.hdr'}',sprintf('Find CT images for %s',template.subject));

    

if nargin < 3
    sdthresh = []; %level at which to threshold z-transformed image
end

if nargin < 4
    sm = 1.2; %smoothing kernel in mm
end

[~,fn,ext] = fileparts(fname);
fname =fullfile(path,fname);

%Load image data
switch lower(ext)
    case {'.hdr'}
        [IM,hdat] = loadimg(fname);

    case {'.nii'}
        [IM,hdat] = readnifti(fname);    

    otherwise
    error('Unrecognized file extension')
end
    
getmsk=1;
if getmsk
    [fname,path] = uigetfile({'*.nii','*.hdr'}',sprintf('Find brain masl for %s',template.subject));
    [~,fn,ext]=fileparts(fname);
    if fname == 0
        MSK = 1;
    else   
        fname = fullfile(path,fname);
        switch lower(ext)
            case {'.hdr'}
                [MSK,mskhdat] = loadimg(fname);

            case {'.nii'}
                [MSK,mskhdat] = readnifti(fname);    

            otherwise
            error('Unrecognized file extension')
        end
    end
else
    MSK=1;
end

IM = IM.*MSK;

voxdim = double([hdat.pixdim1,hdat.pixdim2,hdat.pixdim3]);
imdim = double([hdat.dim1,hdat.dim2,hdat.dim3]);

%%% Compute smoothing kernel
sms = sm./voxdim;

for i = 1:3
    gs{i} = gausswin(ceil(sms(i)*2.5*2));
    if mod(ceil(sms(i)*2.5*2),2) == 1
        gs{i} = [0;gs{i}];
    end
end

fprintf('\nSmoothing...')
G12 = gs{1}*gs{2}';

G = repmat(permute(gs{3},[2 3 1]),size(G12)).*repmat(G12,[1 1 length(gs{3})]);
G = G./sum(G(:));
% fg = fftn(G);

%%% Do quick convolution in fourier space
GG = zeros(size(IM));
GG([1:size(G,1)./2,end-size(G,1)./2+1:end],...
    [1:size(G,2)./2,end-size(G,2)./2+1:end],...
    [1:size(G,3)./2,end-size(G,3)./2+1:end]) = fftshift(G);
FIM = fftn(IM);
FGG = fftn(GG);
SM = ifftn(FIM.*conj(FGG));  % smoothed image
fprintf('Done')

clear GG FIM FGG 


%%% Apply z-score transformation to image
imsc = zeros(size(SM));
imsc(:) = zscore(SM(:));

% ZC = zeroscross(imsc.*(imsc>sdthresh));
ZC = zeroscross(imsc);
[I,J,K] = meshgrid(1:imdim(1),1:imdim(2),1:imdim(3));
IJK = [I(:),J(:),K(:)];
IJK(:,4) = 1;

%transformation matrix from voxels to physical units
trm = double(hdat.vox2unit');
trm(4,4) = 1;


% XYZ = [I(:),J(:),K(:)];


pkijk = IJK(ZC(:),:);
% pkijk(:,end+1) = 1;
pks = pkijk*trm;
pkval = imsc(ZC);

if isempty(sdthresh)
    sdthresh = setthreshold(pks,pkval);
end

pks = pks(pkval>sdthresh,:);



for i = 1:length(template.group)
    template.group(i).XY(:,4) = [0 0 0 1]*hdat.vox2unit'*[ 0 0 1]';
end

T = template.T;
TXYs = {template.group.XY};

% [T,TXYs] = template;

% Use simulate anealing to match peaks to distance template
% tempgroups = {[1 2],[3:length(TXYs)]};

%%

for i = 1:length(template.group)
    Q = TXYs(i);
    Q{1}(:,1) = Q{1}(:,1)- mean(Q{1}(:,1));
    [pksfilt,pk2contacts_locked] = manualmap(cat(1,Q{1}),pks,template.group(i).label); %Peak to contact mapping

    pkD = dfun(pksfilt,pksfilt);

    T = template.group(i).Ds;
%     T(end+1,:) = -1; T(:,end+1) = -1;
    pkc = simanneal(pkD,T,Q,pksfilt,pk2contacts_locked); %Peak to contact mapping
    
%     pk2contacts_locked=pk2contacts_locked|pkc>0;
    
    [srt,srti] = sort(pkc);
    pkmaps(i).pkc = srt(srt>0);
    pkmaps(i).pkssrt = pks(srti(srt>0),:);
   
end

%%

XYZ = IJK*trm;
decim = 4;
X = reshape(XYZ(:,1),size(I));
x = squeeze(X(1,1:decim:end,1));
Y = reshape(XYZ(:,2),size(I));
y = squeeze(Y(1:decim:end,1,1));
Z = reshape(XYZ(:,3),size(I));
z = squeeze(Z(1,1,1:decim:end));
Zdecim = Z(1:decim:end,1:decim:end,1:decim:end);
Xdecim = X(1:decim:end,1:decim:end,1:decim:end);
Ydecim = Y(1:decim:end,1:decim:end,1:decim:end);
IMdecim = IM(1:decim:end,1:decim:end,1:decim:end);
imnrm = (IMdecim-min(IMdecim(:)))./range(IMdecim(:));

figure
 plot3(pksfilt(:,1),pksfilt(:,2),pksfilt(:,3),'ro')

 hold on,
slz = 100/decim;
imh(1) =surf(Xdecim(:,:,slz),Ydecim(:,:,slz),Zdecim(:,:,slz),imnrm(:,:,slz));
hold on
sly = 200/decim;
imh(2) =surf(Xdecim(:,sly,:),Ydecim(:,sly,:),squeeze(Zdecim(:,sly,:)),squeeze(imnrm(:,sly,:)));
slx = 200/decim;
imh(3) =surf(squeeze(Xdecim(slx,:,:)),squeeze(Ydecim(slx,:,:)),squeeze(Zdecim(slx,:,:)),squeeze(imnrm(slx,:,:)));
shading flat
colormap gray
set(imh,'facealpha',.8)
