function in = isPointInTriangle(point, tri_v1, tri_v2, tri_v3)
% ISPOINTINTRIANGLE Checks if a point is inside a triangle using barycentric coordinates.
%
% Inputs:
%   point  - The point to check (1x3 vector).
%   tri_v1, tri_v2, tri_v3 - The three vertices of the triangle (1x3 vectors).
%
% Output:
%   in - True if the point is inside the triangle, false otherwise.

    v0 = tri_v2 - tri_v1;
    v1 = tri_v3 - tri_v1;
    v2 = point - tri_v1;

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    in = (u >= 0) && (v >= 0) && (u + v < 1);
end