clear all; clc

image_file = 'H:\Desktop\WLTQT1_correctOrientation.nii';

V = niftiread(image_file);
info = niftiinfo(image_file);

xyz_loc = [126 98 60];
% ML, AP (larger value = more anterior), DV

x = xyz_loc(1);
y = xyz_loc(2);
z = xyz_loc(3);

figure('Renderer','painters','Position',[200 200 1500 350]); hold on;

subplot(1,3,1)  % Sagittal
imagesc(rot90(squeeze(V(70,:,:))))
set(gca,'YDir','reverse')
title('Sagittal')
axis image off
vline(y,'r-')

subplot(1,3,2); hold on  % Coronal
imagesc(rot90(squeeze(V(:,y,:))))
set(gca,'YDir','reverse')
title('Coronal')
axis image off
scatter(info.ImageSize(1)-x, info.ImageSize(3)-z, 20, 'r', 'filled')

subplot(1,3,3); hold on  % Axial
imagesc(rot90(squeeze(V(:,:,z))))
scatter(info.ImageSize(1)-x, info.ImageSize(2)-y, 20, 'r', 'filled')
set(gca,'YDir','reverse')
title('Axial')
axis image off

colormap gray
