clear all; clc; close all;

[ephysLog, stimulusLog, spike_log] = import_exp_map();
ephysLog = ephysLog(strcmp(ephysLog.monkey,'walt'),:);

% Coordinates are already in voxel space
vox_x = round(ephysLog.nan_x);           % ML
vox_y = round(ephysLog.nan_y);           % AP
vox_z = 167-3-round(ephysLog.nan_z_target); % DV (flipped)
area = ephysLog.area_label_sec;
frontal_mask  = strcmp(area, 'frontal');
auditory_mask = strcmp(area, 'auditory');

%% Load MRI
nii = load_nii('WLTQT1_correctOrientation.nii');
img = double(nii.img);
img_flipped = flip(img, 1);
img_flipped = flip(img_flipped, 3);

dims = size(img_flipped);

% Clamp to image dimensions
vox_x = max(1, min(vox_x, dims(1)));
vox_y = max(1, min(vox_y, dims(2)));
vox_z = max(1, min(vox_z, dims(3)));

%% Slice indices
frontal_y_slice  = round(mode(vox_y(frontal_mask)));
auditory_y_slice = round(mode(vox_y(auditory_mask)));

% Zoom windows
zoom_pad = 12;
clamp_range = @(c, pad, mx) max(1, c-pad) : min(mx, c+pad);

frontal_xr  = clamp_range(round(median(vox_x(frontal_mask))),  zoom_pad, dims(1));
frontal_zr  = clamp_range(round(median(vox_z(frontal_mask))),  zoom_pad, dims(3));
auditory_xr = clamp_range(round(median(vox_x(auditory_mask))), zoom_pad, dims(1));
auditory_zr = clamp_range(round(median(vox_z(auditory_mask))), zoom_pad, dims(3));

tol = 5; % voxels either side of display slice

%% Shared zoom ranges for sagittal (slice along X, plot Y vs Z)
frontal_yr  = clamp_range(round(median(vox_y(frontal_mask))),  zoom_pad, dims(2));
frontal_zr2 = clamp_range(round(median(vox_z(frontal_mask))),  zoom_pad, dims(3));
auditory_yr  = clamp_range(round(median(vox_y(auditory_mask))), zoom_pad, dims(2));
auditory_zr2 = clamp_range(round(median(vox_z(auditory_mask))), zoom_pad, dims(3));

frontal_x_slice  = round(mode(vox_x(frontal_mask)));
auditory_x_slice = round(mode(vox_x(auditory_mask)));

%% Figure 1: Coronal
figure('Color','k','Position',[100 100 900 900]);

% ── Subplot 1: Frontal full coronal ───────────────────────────────────
ax1 = subplot(2,2,1);
imagesc(squeeze(img_flipped(:, frontal_y_slice, :))');
colormap(ax1, gray); axis image off; hold on;
in_slice = abs(vox_y - frontal_y_slice) <= tol & frontal_mask;
scatter(vox_x(in_slice), vox_z(in_slice), 12, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
title('Frontal', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 2: Auditory full coronal ──────────────────────────────────
ax2 = subplot(2,2,2);
imagesc(squeeze(img_flipped(:, auditory_y_slice, :))');
colormap(ax2, gray); axis image off; hold on;
in_slice = abs(vox_y - auditory_y_slice) <= tol & auditory_mask;
scatter(vox_x(in_slice), vox_z(in_slice), 12, 'b', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
title('Auditory', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 3: Frontal zoomed coronal ─────────────────────────────────
ax3 = subplot(2,2,3);
imagesc(squeeze(img_flipped(frontal_xr, frontal_y_slice, frontal_zr))');
colormap(ax3, gray); axis image off; hold on;
in_slice = abs(vox_y - frontal_y_slice) <= tol & frontal_mask;
scatter(vox_x(in_slice) - frontal_xr(1) + 1, ...
        vox_z(in_slice) - frontal_zr(1) + 1, ...
        20, 'r', 'MarkerEdgeColor', 'w');
title('Frontal (zoom)', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 4: Auditory zoomed coronal ────────────────────────────────
ax4 = subplot(2,2,4);
imagesc(squeeze(img_flipped(auditory_xr, auditory_y_slice, auditory_zr))');
colormap(ax4, gray); axis image off; hold on;
in_slice = abs(vox_y - auditory_y_slice) <= tol & auditory_mask;
scatter(vox_x(in_slice) - auditory_xr(1) + 1, ...
        vox_z(in_slice) - auditory_zr(1) + 1, ...
        20, 'b', 'MarkerEdgeColor', 'w');
title('Auditory (zoom)', 'Color', 'w');
set(gca, 'YDir', 'reverse');

%% Figure 2: Sagittal
figure('Color','k','Position',[100 100 900 900]);

% ── Subplot 1: Frontal full sagittal ──────────────────────────────────
ax5 = subplot(2,2,1);
imagesc(squeeze(img_flipped(frontal_x_slice, :, :))');
colormap(ax5, gray); axis image off; hold on;
in_slice = abs(vox_x - frontal_x_slice) <= tol & frontal_mask;
scatter(vox_y(in_slice), vox_z(in_slice), 12, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
title('Frontal', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 2: Auditory full sagittal ─────────────────────────────────
ax6 = subplot(2,2,2);
imagesc(squeeze(img_flipped(auditory_x_slice, :, :))');
colormap(ax6, gray); axis image off; hold on;
in_slice = abs(vox_x - auditory_x_slice) <= tol & auditory_mask;
scatter(vox_y(in_slice), vox_z(in_slice), 12, 'b', 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
title('Auditory', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 3: Frontal zoomed sagittal ────────────────────────────────
ax7 = subplot(2,2,3);
zoom_img = squeeze(img_flipped(frontal_x_slice, frontal_yr, frontal_zr2))';
imagesc(zoom_img);
colormap(ax7, gray); axis image off; hold on;
in_slice = abs(vox_x - frontal_x_slice) <= tol & frontal_mask;
scatter(vox_y(in_slice) - frontal_yr(1) + 1, ...
        vox_z(in_slice) - frontal_zr2(1) + 1, ...
        20, 'r', 'MarkerEdgeColor', 'w');
title('Frontal (zoom)', 'Color', 'w');
set(gca, 'YDir', 'reverse');

% ── Subplot 4: Auditory zoomed sagittal ───────────────────────────────
ax8 = subplot(2,2,4);
zoom_img = squeeze(img_flipped(auditory_x_slice, auditory_yr, auditory_zr2))';
imagesc(zoom_img);
colormap(ax8, gray); axis image off; hold on;
in_slice = abs(vox_x - auditory_x_slice) <= tol & auditory_mask;
scatter(vox_y(in_slice) - auditory_yr(1) + 1, ...
        vox_z(in_slice) - auditory_zr2(1) + 1, ...
        20, 'b', 'MarkerEdgeColor', 'w');
title('Auditory (zoom)', 'Color', 'w');
set(gca, 'YDir', 'reverse');