function are_overlapping = check_scatterers_overlap(center1, dim1, center2, dim2)
% check_scatterers_overlap: Checks if the 3D bounding boxes of two scatterers overlap.
%
% Inputs:
%   center1, center2: [x, y, z] center coordinates of the scatterers.
%   dim1, dim2:       [dx, dy, dz] dimensions (depth, width, height) of the scatterers.
%
% Output:
%   are_overlapping:  A boolean (true if they overlap, false otherwise).

    % Calculate the min and max coordinates for the first scatterer
    min1 = center1 - dim1 / 2;
    max1 = center1 + dim1 / 2;
    
    % Calculate the min and max coordinates for the second scatterer
    min2 = center2 - dim2 / 2;
    max2 = center2 + dim2 / 2;
    
    % Check for overlap on each of the three axes
    x_overlap = (max1(1) > min2(1)) && (max2(1) > min1(1));
    y_overlap = (max1(2) > min2(2)) && (max2(2) > min1(2));
    z_overlap = (max1(3) > min2(3)) && (max2(3) > min1(3));
    
    % Overlap exists only if there is overlap on all three axes simultaneously
    are_overlapping = x_overlap && y_overlap && z_overlap;
end