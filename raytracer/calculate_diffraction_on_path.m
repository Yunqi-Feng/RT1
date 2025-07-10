function [total_diffraction_loss, diffraction_points] = calculate_diffraction_on_path(start_point, end_point, edges, freq)
% calculate_diffraction_on_path: Checks for diffraction on a path segment.
%
% Inputs:
%   start_point: Coordinates of the start of the path segment.
%   end_point:   Coordinates of the end of the path segment.
%   edges:       A list of all edges in the environment.
%   lambda:      Wavelength of the signal.
%
% Outputs:
%   total_diffraction_loss: The total diffraction loss for the segment in dB.
%   diffraction_points:     The coordinates of the points of diffraction.
    lambda = 3e8 / freq;
    total_diffraction_loss = 0;
    diffraction_points = [];
    path_vector = end_point - start_point;
    path_length = norm(path_vector);

    for i = 1:size(edges, 1)
        edge = [edges(i, 1:3); edges(i, 4:6)];
        
        % Find the point on the edge closest to the path
        p1 = start_point;
        p2 = end_point;
        e1 = edge(1, :);
        e2 = edge(2, :);
        
        % Vector perpendicular to both the path and the edge
        cross_prod = cross(path_vector, e2 - e1);
        
        if norm(cross_prod) > 1e-6 % If the path and edge are not parallel
            t = dot(cross(e1 - p1, e2 - e1), cross_prod) / norm(cross_prod)^2;
            s = dot(cross(e1 - p1, path_vector), cross_prod) / norm(cross_prod)^2;
            
            if t >= 0 && t <= 1 && s >= 0 && s <= 1
                closest_point_on_path = p1 + t * path_vector;
                closest_point_on_edge = e1 + s * (e2 - e1);
                
                h = norm(closest_point_on_edge - closest_point_on_path);
                d1 = norm(closest_point_on_path - start_point);
                d2 = path_length - d1;
                
                % Check if the edge is in the first Fresnel zone
                fresnel_radius = sqrt(lambda * d1 * d2 / (d1 + d2));
                
                if h > 0 && h < fresnel_radius
                    loss = calculate_ked_loss(h, d1, d2, lambda);
                    if loss > total_diffraction_loss % Consider the worst-case diffraction
                        total_diffraction_loss = loss;
                        diffraction_points = closest_point_on_edge;
                    end
                end
            end
        end
    end
end