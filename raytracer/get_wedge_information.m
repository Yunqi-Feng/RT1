function wedges = get_wedge_information(valid_edges, CADop)
% get_wedge_information: Finds the wedge information (two planes, materials) for a pre-filtered list of valid edges.
% This function is designed to work with the output of generate_edges.m.
%
% Inputs:
%   valid_edges: An Mx6 matrix of edge coordinates, where M is the number of valid edges
%                (This is the output of generate_edges.m).
%   CADop:       The full, original CAD output from xmlreader, containing all triangles and their properties.
%
% Output:
%   wedges:      A struct array where each element corresponds to a valid edge and contains:
%                .edge_v1, .edge_v2: The two vertices of the shared edge.
%                .plane1, .plane2:   The plane equations of the two faces.
%                .mat_idx1, .mat_idx2: The material library indices for each plane.

    % Initialize the output structure
    wedges = struct('edge_v1', {}, 'edge_v2', {}, 'plane1', {}, 'plane2', {}, 'mat_idx1', {}, 'mat_idx2', {}, 'n', {});
    
    % Use a tolerance for comparing floating-point vertex coordinates
    tolerance = 1e-5;

    % --- Main Loop to Find Wedges for Each Valid Edge ---
    for i = 1:size(valid_edges, 1)
        current_edge_v1 = valid_edges(i, 1:3);
        current_edge_v2 = valid_edges(i, 4:6);
        
        found_planes_indices = [];
        
        % Search through the entire CADop to find triangles that contain this edge
        for j = 1:size(CADop, 1)
            triangle_verts = [CADop(j, 1:3); CADop(j, 4:6); CADop(j, 7:9)];
            
            % Check if both vertices of the current edge are part of this triangle
            has_v1 = any(all(abs(triangle_verts - current_edge_v1) < tolerance, 2));
            has_v2 = any(all(abs(triangle_verts - current_edge_v2) < tolerance, 2));
            
            if has_v1 && has_v2
                % This triangle (j) is one face of the wedge
                found_planes_indices(end+1) = j;
            end
        end
        
        % A valid wedge must be formed by exactly two planes.
        % If we found two planes, store the wedge information.
        if length(found_planes_indices) == 2
            idx1 = found_planes_indices(1);
            idx2 = found_planes_indices(2);
            
            new_wedge.edge_v1 = current_edge_v1;
            new_wedge.edge_v2 = current_edge_v2;
            

            plane1_eq = CADop(idx1, 10:13);
            plane2_eq = CADop(idx2, 10:13);


            new_wedge.plane1 = CADop(idx1, 10:13);
            new_wedge.plane2 = CADop(idx2, 10:13);
            
            new_wedge.mat_idx1 = CADop(idx1, 14);
            new_wedge.mat_idx2 = CADop(idx2, 14);
            

            % Calculate wedge angle
            normal1 = plane1_eq(1:3);
            normal2 = plane2_eq(1:3);
            wedge_angle = acos(dot(normal1, normal2));
            
            % Calculate wedge parameter 'n'
            new_wedge.n = (2*pi - wedge_angle) / pi;

            wedges(end+1) = new_wedge;
        else
            % This case might occur for edges on the boundary of an object
            % that don't form a closed wedge with another surface.
            %fprintf('Note: Edge %d does not form a clean 2-plane wedge. Found %d planes. Skipping.\n', i, length(found_planes_indices));
        end
    end
    
    %fprintf('Found wedge information for %d valid edges.\n', length(wedges));
end