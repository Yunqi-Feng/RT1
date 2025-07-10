function [intersection_point, check] = plane_line_intersect(plane_normal, point_on_plane, line_point, line_direction)
% PLANE_LINE_INTERSECT Calculates the intersection of a line and a plane in 3D.
%
% This function finds the point where a line, defined by a point and a
% direction vector, intersects with a plane, defined by its normal vector
% and a point on the plane.
%
% Inputs:
%   plane_normal    - A 1x3 vector representing the normal to the plane.
%   point_on_plane  - A 1x3 vector representing a known point on the plane.
%   line_point      - A 1x3 vector representing a known point on the line.
%   line_direction  - A 1x3 vector representing the direction of the line.
%
% Outputs:
%   intersection_point - The 1x3 coordinates of the intersection point.
%                        Returns an empty array if there is no unique intersection.
%   check              - An integer flag:
%                        - 1: A unique intersection point was found.
%                        - 0: The line is parallel to the plane (no intersection).
%                        - 2: The line lies within the plane (infinite intersections).

    % Ensure inputs are row vectors
    plane_normal = plane_normal(:)';
    point_on_plane = point_on_plane(:)';
    line_point = line_point(:)';
    line_direction = line_direction(:)';

    ndotv = dot(plane_normal, line_direction);

    % Check if the line is parallel to the plane
    if abs(ndotv) < 1e-6
        % Check if the line lies within the plane
        if abs(dot(plane_normal, line_point - point_on_plane)) < 1e-6
            check = 2; % Line is in the plane
            intersection_point = []; % No unique intersection
        else
            check = 0; % Line is parallel and not in the plane
            intersection_point = []; % No intersection
        end
        return;
    end

    % Find the intersection point using the formula:
    % t = dot(plane_normal, point_on_plane - line_point) / dot(plane_normal, line_direction)
    % Intersection_Point = line_point + t * line_direction
    
    t = dot(plane_normal, point_on_plane - line_point) / ndotv;
    
    intersection_point = line_point + t * line_direction;
    check = 1; % Unique intersection found

end