function dist = euclidean_distance_3d(p1, p2)
%EUCLIDEAN_DISTANCE_3D Computes the Euclidean distance between two 3D points
%
%   dist = euclidean_distance_3d(p1, p2)
%
%   Inputs:
%       p1 - 1x3 vector representing the first 3D point [x1, y1, z1]
%       p2 - 1x3 vector representing the second 3D point [x2, y2, z2]
%
%   Output:
%       dist - Euclidean distance between the two points

    % Ensure the inputs are 1x3 vectors
    assert(all(size(p1) == [1 3]), 'p1 must be a 1x3 vector');
    assert(all(size(p2) == [1 3]), 'p2 must be a 1x3 vector');

    % Compute Euclidean distance
    dist = sqrt(sum((p1 - p2).^2));
end
