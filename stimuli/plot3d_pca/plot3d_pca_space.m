clear all; clc

load('pc_data_in.mat')

outfile_name = 'frontal.mp4';

% Extract the first three principal components from the frontal dataset
x = pc_data_in(:,1);   % PC1 values
y = pc_data_in(:,2);   % PC2 values
z = pc_data_in(:,3);   % PC3 values (for 3D plotting)

% Define timestamps (in ms) of sound events
event_times = [0, 563, 1126, 1689, 2252];

% Convert sound times to frame indices by adding an offset (e.g., 100ms baseline)
event_times_idx = event_times + 100;

% Define indices at which the plot will be updated (every 10 steps)
plot_idx = 2:10:length(x);

% Axis range offset for projection lines
axis_limits = 5;

% Find the closest indices in plot_idx corresponding to event times
for event_i = 1:length(event_times)
    [~, plot_sound_idx(event_i)] = min(abs(plot_idx - event_times_idx(event_i)));
end

% Create a video writer object for output video (MPEG-4 format)
v = VideoWriter(outfile_name, 'MPEG-4');
v.FrameRate = 30;  % Set video frame rate to 30 fps
open(v);  % Open the video file for writing

% Create and configure the figure
figuren('Renderer', 'painters', 'Position', [200 200 700 700]); hold on;
subplot(3,1,[1 2 3]); hold on  % Use most of the figure space for a single 3D plot
view(30, 28);  % Set 3D view angle
grid on
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

% Initialize graphic object handles for dynamic updating
pc1_line = []; pc2_line = []; pc3_line = []; pcall_line = [];

% Loop over time steps and plot data incrementally
for k = 2:10:length(x)
    % Delete previous lines to avoid overplotting
    if isgraphics(pc1_line); delete(pc1_line); end
    if isgraphics(pc2_line); delete(pc2_line); end
    if isgraphics(pc3_line); delete(pc3_line); end
    if isgraphics(pcall_line); delete(pcall_line); end

    % Plot the full colored 3D trajectory up to time step k
    pcall_line = color_line3(x(1:k), y(1:k), z(1:k), 1:k);

    % Function: 

    % Plot projections of trajectory onto PC2-PC3, PC1-PC3, PC1-PC2 planes
    pc1_line = plot3(ones(1,k)*-axis_limits, y(1:k), z(1:k), 'Color', [1 1 1]*0.8);
    pc2_line = plot3(x(1:k), ones(1,k)*axis_limits, z(1:k), 'Color', [1 1 1]*0.8);
    pc3_line = plot3(x(1:k), y(1:k), ones(1,k)*-axis_limits, 'Color', [1 1 1]*0.8);

    % Highlight points where a sound occurred
    if ismember(k, plot_idx(plot_sound_idx))
        scatter3(x(k), y(k), z(k), ...
            100, [0 0 0], '^', 'filled');  % Black triangle marker
    end

    % Set plot title with current time (ms)
    title([int2str(k-100) 'ms'], 'Units', 'normalized', 'Position', [0.1, -0.05, 0]);

    % Set axis limits
    axis([-1 1 -1 1 -1 1] * axis_limits);

    % Capture the current frame and write it to the video
    frame = getframe(gcf);
    writeVideo(v, frame);

    % Update figure
    drawnow;
    pause(0.01);  % Small pause for visualization
end

% Rotate the final 3D plot around the z-axis for a 360° view
for azi = 1:360
    view(30 + azi, 28);  % Increment azimuthal angle
    pause(0.01);         % Pause for visual smoothness
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Finalize and close the video file
close(v);
