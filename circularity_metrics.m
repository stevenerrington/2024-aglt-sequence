function C = circularity_metrics(PC1, PC2)
% CIRCULARITY_METRICS Computes common circularity measures for 2D trajectories.
%
% Inputs:
%   PC1, PC2  = trajectory coordinates (vectors of same length)
%
% Outputs (struct C):
%   C.R        = mean resultant length (1 = perfect circular)
%   C.V        = circular variance       (0 = perfect circular)
%   C.CV_r     = coefficient of variation of radius
%   C.ellipse  = ellipticity circularity index (1=circle, 0=line)
%

% --- 1. Angle-based circularity (mean resultant length) ---
theta = atan2(PC2, PC1);
C.R = abs(mean(exp(1i * theta)));     % mean resultant length
C.V = 1 - C.R;                         % circular variance

% --- 2. Radius stability ---
r = sqrt(PC1.^2 + PC2.^2);
C.CV_r = std(r) / mean(r);             % low = circular

% --- 3. Ellipticity (ratio of eigenvalues of covariance) ---
X = [PC1(:) PC2(:)];
X = X - mean(X);                       % center trajectory
COV = cov(X);
[evecs, evals] = eig(COV);
evals = sort(diag(evals), 'descend');  % λ1 ≥ λ2

lambda1 = evals(1);
lambda2 = evals(2);

% Circularity index: 1 = perfect circle, 0 = fully elongated
C.ellipse = 1 - (lambda1 - lambda2) / (lambda1 + lambda2);
end
