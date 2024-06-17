

% ----------- SVN REVISION INFO ------------------
% $URL: https://saccade.neurosurgery.uiowa.edu/svn/KovachToolbox/devel/SubNuclear/do_registration.m $
% $Revision: 401 $
% $Date: 2013-10-28 10:33:55 -0500 (Mon, 28 Oct 2013) $
% $Author: ckovach $
% ------------------------------------------------

%% Get Subject ID
sid = inputdlg('Which subject?');
sid= sid{1};

ddir = sprintf('%sfsl.anat',sid);

%% Subcortical parcellation
% Calls FSL to 
run_fsl_script

%% Load Data

KTcheck
postop = volumeview(sprintf('%s PostOp',sid),'post_op_aligned.nii.gz',ddir);
preop = volumeview(sprintf('%s PreOp',sid),'T1.nii.gz',ddir);
postop.intensity_range = [0 400];
mnih = readnifti('MNI152_T1_1mm.nii',true);

preop.sisters= postop;

%%% Compute the voxel to MNI coordinate transformation
T = textread(fullfile(ddir,'T1_to_MNI_lin.mat'))'; %#ok<REMFF1>
tr1 = transforms('trmat',T(1:4,:),'label','T12MNI');
tr2 = transforms('trmat',double(mnih.vox2unit)','label','MNI2mm');
tr = tr1*tr2;
tr.label = 'vox2MNImm';
postop.addtransform(tr);
preop.addtransform(tr);







helpdlg('Check the alignment of the images and select matching control points on the preop and postop brains')
%% Atlas Coregistration for Amygdala
% 
% After you're satisfied with the amygdala parcellation, meshes for amygdalae are loaded from the file,  mai_template_mni_aligned.mat, and subnuclei are loaded from maiwarp2mni.mat
%  These are extracted from “Atlas of the Human Brain” (Mai 2008). The variables, newmeshR and newmeshL contain data for right and left amygdalae, respectively, as TriRep objects  (see matlab help for details on the TriRep class). The vertices in newmeshL and newmeshR are matched to those in the FSL template mesh for the amygdala. Coregistration uses thin plate spline warping to align the vertices in the FSL generated mesh with the atlas-derived mesh. The necessary steps are carried out in do_registration. The resulting warping function is then applied to the vertices of meshes for each subnucleus in order to transform the subnuclei into the subject's image space. 
% 

structs = {'L Amyg' , fullfile(ddir,'first_results/T1_first-L_Amyg_first.vtk')
           'R Amyg' , fullfile(ddir,'first_results/T1_first-R_Amyg_first.vtk')};

for i = 1:length(structs)
%     postop.addmesh(structs{i,2},structs{i,1})
     preop.addmesh(structs{i,2},structs{i,1})
end

helpdlg('Now adjust the amygdalae meshes')


%%
% helpdlg('Now adjust the meshes')
% 
% for i = 3:length(preop.meshes)
%    postop.addmesh(preop.meshes(i));
%    
% end


%% Save progress
save(fullfile(ddir,sprintf('%s_volview',sid)),'preop','postop'),
%% Next apply TPS warping


helpdlg('Select control points on each image. Make sure they are in the same order')
x1 = cat(1,preop.points.coord); x1(:,4) = 1;
x2 = cat(1,postop.points.coord); x2(:,4) = 1;



mwarp = tpswarp(x1(:,1:3),x2(:,1:3)); % this warps from preop to postop
imwarp = tpswarp(x2(:,1:3),x1(:,1:3)); % this from postop to preop
% mx2tps  = mwarp(mx1(:,1:3));

XL = preop.meshes(end-1).trirep.X; XL(:,4)=1;
% mlinL = TriRep(preop.meshes(1).trirep.Triangulation,XL*TlinL(:,1:3));
mtpsL = TriRep(preop.meshes(end-2).trirep.Triangulation,mwarp(XL(:,1:3)));
XR = preop.meshes(end).trirep.X; XR(:,4)=1;
% mlinR = TriRep(preop.meshes(2).trirep.Triangulation,XR*TlinR(:,1:3));
mtpsR = TriRep(preop.meshes(end-1).trirep.Triangulation,mwarp(XR(:,1:3)));


%%
figure(postop.fig)
helpdlg('Now add the contact locations')
%% Project contacts to preop brain
for i =length(preop.points)+1:length(postop.points);
    pt= postop.points(i);
    
    preop.addpoint(pt.label,imwarp(pt.coord));
end

%%  Project the mesh from the preop brain onto the postop as an additional sanity check
% postop.addmesh(mlinL,'LA Lin');
postop.addmesh(mtpsL,'LA TPS');
% postop.addmesh(mlinR,'RA Lin');
postop.addmesh(mtpsR,'RA TPS');

%% Now load the mai atlas subnuclei
load maiwarp2mni  % Atlas warped into MNI space based on FSL meshes
load mai_template_mni_aligned newmeshR newmeshL tmpl2maimatR tmpl2maimatL % Atlas with vertices matcehd to FSL vertices

%%
tr = preop.transforms(3); %This sould be the vox to MNI transform


warpfunL  = tpswarp(tr.itr(newmeshL.X),preop.meshes(3).trirep.X,.1);  % Left amygdala

warpfunR = tpswarp(tr.itr(newmeshR.X),preop.meshes(4).trirep.X,.1); % Right amygdala

[iwarpfunL,AL]   = tpswarp(tr.tr(preop.meshes(3).trirep.X),newmeshL.X*10);
[iwarpfunR,AR] = tpswarp(tr.tr(preop.meshes(4).trirep.X),newmeshR.X*10);

%%% Affine transformation into template space (so we can match the axial
%%% view to the nearest page in the template).
AL(4,4) = 1; AR(4,4) = 1;
tmpl2maimatL(4,4) = 1; tmpl2maimatR(4,4) = 1;
preop.addtransform(tr.trmat*AL*tmpl2maimatL,'Left Amyg. vox 2 atlas');
preop.addtransform(tr.trmat*AR*tmpl2maimatR,'Right Amyg. vox 2 atlas');


kp = find(arrayfun(@(x)~isempty(x.warped),maiwarpL));
kp(end) = [];

colors = {'k','r','g','y','m','c','c','k','b','b'};
for i = 1:length(kp)   
    wrpR = TriRep(maiwarpR(kp(i)).warped.Triangulation,warpfunR(tr.itr(maiwarpR(kp(i)).warped.X)));
    preop.addmesh(wrpR,sprintf('R %s',maiwarpR(kp(i)).label)); 
    
end
for i = 1:length(kp)   
    wrpL = TriRep(maiwarpL(kp(i)).warped.Triangulation,warpfunL(tr.itr(maiwarpL(kp(i)).warped.X)));
    preop.addmesh(wrpL,sprintf('L %s',maiwarpL(kp(i)).label));
end

