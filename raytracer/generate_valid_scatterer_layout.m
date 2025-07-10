function [x_cen, y_cen, z_cen, x_dep, y_dep, z_dep] = generate_valid_scatterer_layout(desired_scatterer_num, confg_dim, scatterer_dim, scatterer_dis, NodeLocTx, NodeLocRx, IndoorSwitch, scatterer_confg, scatterer_height_max)
% generate_valid_scatterer_layout: Places a desired number of non-overlapping scatterers.

    %fprintf('--> Starting iterative placement for %d scatterers...\n', desired_scatterer_num);
    max_attempts_per_scatterer = 200; % Max tries to place a single scatterer

    % These arrays will store the final, validated scatterer properties
    x_cen_final = []; y_cen_final = []; z_cen_final = [];
    x_dep_final = []; y_dep_final = []; z_dep_final = [];

    for i = 1:desired_scatterer_num
        is_placed_successfully = false;
        attempt_count = 0;

        while ~is_placed_successfully && attempt_count < max_attempts_per_scatterer
            % 1. Generate ONE candidate scatterer's properties
            [cand_xc, cand_yc, cand_zc, cand_xd, cand_yd, cand_zd, ~] = scatterer_config(confg_dim, 1, scatterer_dim, scatterer_dis, NodeLocTx, NodeLocRx, IndoorSwitch, [], scatterer_height_max);
            
            % 2. Check for overlap with already placed scatterers
            has_scatterer_overlap = false;
            for j = 1:length(x_cen_final)
                if check_scatterers_overlap([cand_xc, cand_yc, cand_zc], [cand_xd, cand_yd, cand_zd], [x_cen_final(j), y_cen_final(j), z_cen_final(j)], [x_dep_final(j), y_dep_final(j), z_dep_final(j)])
                    has_scatterer_overlap = true;
                    break;
                end
            end
            
            % 3. Check for overlap with Tx/Rx nodes
            tx_is_inside = (NodeLocTx(1) >= cand_xc - 0.5*cand_xd && NodeLocTx(1) <= cand_xc + 0.5*cand_xd && ...
                            NodeLocTx(2) >= cand_yc - 0.5*cand_yd && NodeLocTx(2) <= cand_yc + 0.5*cand_yd && ...
                            NodeLocTx(3) >= cand_zc - 0.5*cand_zd && NodeLocTx(3) <= cand_zc + 0.5*cand_zd);
            rx_is_inside = (NodeLocRx(1) >= cand_xc - 0.5*cand_xd && NodeLocRx(1) <= cand_xc + 0.5*cand_xd && ...
                            NodeLocRx(2) >= cand_yc - 0.5*cand_yd && NodeLocRx(2) <= cand_yc + 0.5*cand_yd && ...
                            NodeLocRx(3) >= cand_zc - 0.5*cand_zd && NodeLocRx(3) <= cand_zc + 0.5*cand_zd);

            % 4. If no overlaps are found, the placement is valid
            if ~has_scatterer_overlap && ~tx_is_inside && ~rx_is_inside
                % Success! Add the candidate to our final list.
                x_cen_final(end+1) = cand_xc; y_cen_final(end+1) = cand_yc; z_cen_final(end+1) = cand_zc;
                x_dep_final(end+1) = cand_xd; y_dep_final(end+1) = cand_yd; z_dep_final(end+1) = cand_zd;
                is_placed_successfully = true;
                %fprintf('    Placed scatterer %d/%d successfully.\n', i, desired_scatterer_num);
            else
                % Failed placement, try again.
                attempt_count = attempt_count + 1;
            end
        end

        % If we exit the while loop due to max attempts, the scene is too crowded.
        if ~is_placed_successfully
            error('Failed to place scatterer %d after %d attempts. The scene is too crowded. Please reduce the number or size of scatterers.', i, max_attempts_per_scatterer);
        end
    end
    
    % Return the final, validated lists
    x_cen = x_cen_final; y_cen = y_cen_final; z_cen = z_cen_final;
    x_dep = x_dep_final; y_dep = y_dep_final; z_dep = z_dep_final;
end