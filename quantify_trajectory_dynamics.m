%% Trajectory Dynamics Quantification
% X is an n x 3 matrix (positions in state space: [x y z])
% t is an n x 1 time vector (optional, assume dt=1 if not given)

function metrics = quantify_trajectory_dynamics(X, t)

if nargin < 2
    t = (0:size(X,1)-1)'; % default: unit time steps
end

% Ensure column vector for t
t = t(:);

% Time step (assume uniform sampling)
dt = mean(diff(t));

%% 1. Compute Velocity, Acceleration, Jerk
V = gradient(X, dt);           % Velocity (n x 3)
A = gradient(V, dt);           % Acceleration (n x 3)
J = gradient(A, dt);           % Jerk (n x 3)

speed = vecnorm(V,2,2);        % Speed (norm of velocity)
acc_mag = vecnorm(A,2,2);      % Acceleration magnitude
jerk_mag = vecnorm(J,2,2);     % Jerk magnitude

%% 2. Path Length
ds = sqrt(sum(diff(X).^2, 2)); % Incremental path lengths
path_length = sum(ds);

%% 3. Curvature (for 3D trajectories)
% κ = |v x a| / |v|^3
cross_va = cross(V, A, 2);
curvature = vecnorm(cross_va, 2, 2) ./ (speed.^3 + eps);

%% 4. RMS values
v_rms = sqrt(mean(speed.^2));
a_rms = sqrt(mean(acc_mag.^2));
j_rms = sqrt(mean(jerk_mag.^2));

%% 5. Pack results
metrics = struct();
metrics.path_length = path_length;
metrics.speed = speed;
metrics.acceleration = acc_mag;
metrics.jerk = jerk_mag;
metrics.curvature = curvature;
metrics.v_rms = v_rms;
metrics.a_rms = a_rms;
metrics.j_rms = j_rms;
metrics.V = V;
metrics.A = A;
metrics.J = J;


%% Geometric Features
net_displacement = norm(X(end,:) - X(1,:));
efficiency = net_displacement / (path_length + eps);
tortuosity = path_length / (net_displacement + eps);

% Convex hull volume
try
    [~, hull_volume] = convhull(X(:,1), X(:,2), X(:,3));
catch
    hull_volume = NaN; % In case trajectory is collinear
end

% Radius of gyration
centroid = mean(X,1);
R_g = sqrt(mean(sum((X - centroid).^2,2)));

% Principal component analysis
[~,~,latent] = pca(X);
pca_var_ratio = latent / sum(latent);

% Store
metrics.net_displacement = net_displacement;
metrics.efficiency = efficiency;
metrics.tortuosity = tortuosity;
metrics.hull_volume = hull_volume;
metrics.R_g = R_g;
metrics.pca_var_ratio = pca_var_ratio;



end
