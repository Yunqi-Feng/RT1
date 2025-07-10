function [output, multipath, tmp_plane, plane_mat] = ...
    reflectionDiffractionPathGenerator(Tx, Rx, wedges, CADop, vtx, vrx, ...
    frequency, LOS_output,pseudo_LOS_output,polarization, Jones,Ptx, antenna, enablePhase, materialLibrary, ...
    mat_confg, IndoorSwitch, tmp_plane, plane_mat,scattering)
% REFLECTIONDIFFRACTIONPATHGENERATOR Finds paths that undergo one reflection
% then one diffraction.
% ... (rest of the help text)

    output = [];
    multipath = [];
    c = 3e8;
    wavelength = c / frequency;
    indexOutput = 1;

    for i = 1:size(CADop, 1)
        reflecting_plane_eq = CADop(i, 10:13);
        reflecting_plane_verts = [CADop(i, 1:3); CADop(i, 4:6); CADop(i, 7:9)];

        % 1. Create the Image Transmitter
        Tx_image = reflectedImagePointPlane(Tx, reflecting_plane_eq);
        vtx_image = reflectedVelocity(vtx, reflecting_plane_eq);

        % 2. Find Diffraction Paths from the Image Source
        [diffraction_output_from_image, multipath_tmp, tmp_plane, plane_mat,diffraction_loss_dB] = ...
            diffractionPathGenerator_image(Tx_image, Rx, wedges, CADop, ...
            vtx_image, vrx, frequency, LOS_output,pseudo_LOS_output,polarization, Jones,Ptx, antenna, ...
            enablePhase, materialLibrary, mat_confg, IndoorSwitch, ...
            tmp_plane, plane_mat,reflecting_plane_eq);

        % 3. Validate and Reconstruct Each Path Found
        % --- 3. Validate and Reconstruct Each Path Found ---
        if ~isempty(diffraction_output_from_image)
            for j = 1:size(diffraction_output_from_image, 1)
                path_data = diffraction_output_from_image(j, :);
                diffraction_point = multipath_tmp(j,4:6);
        
                % a. Find the reflection point
                % The reflection point is the intersection of the line from the diffraction
                % point to the image source and the reflecting plane.
                dir_vector = Tx_image - diffraction_point;
                [intersection_point, check] = plane_line_intersect(reflecting_plane_eq(1:3), ...
                    reflecting_plane_verts(1,:), diffraction_point, dir_vector);
        
                if check == 0 % No intersection or line is in the plane, which is unlikely
                    continue;
                end
                reflection_point = intersection_point;
        
                % b. Validate the reflection point
                % Check if the reflection point is within the physical triangle.
                is_in_triangle = isPointInTriangle(reflection_point, ...
                    reflecting_plane_verts(1,:), reflecting_plane_verts(2,:), reflecting_plane_verts(3,:));
        
                if ~is_in_triangle
                    continue;
                end
        
                % c. Validate the first leg of the path (Tx -> Reflection Point)
                % The second leg (Reflection Point -> Diffraction -> Rx) is implicitly valid
                % because we found a valid diffraction path from Tx_image.
                is_path_to_reflection_clear = verifyPath(Tx, reflection_point, ...
                    reflection_point - Tx, [0,0,0,0], [0,0,0,0], CADop, 2, false);
        
                if ~is_path_to_reflection_clear
                    continue;
                end
        
                % d. Calculate path properties and losses
                d_tx_refl = norm(reflection_point - Tx);
                d_refl_diff = norm(diffraction_point - reflection_point);
                d_diff_rx = norm(Rx - diffraction_point);
                total_distance = d_tx_refl + d_refl_diff + d_diff_rx;
                delay = total_distance / c;
        
                % Get reflection loss
                % [reflectionLoss_dB, ~] = getDbarcReflectionLoss(Tx, reflection_point, ...
                %     reflecting_plane_eq, CADop(i, 14), materialLibrary, frequency, polarization);
                % [reflectionLoss_dB, ~] = getDbarcReflectionLoss(diffraction_point, reflection_point, ...
                %         reflecting_plane_eq, CADop(i, 14), materialLibrary, frequency, polarization);
                % multipath = [multipath; [Tx, diffraction_point, Rx]];
                multipathtmp = [1,diffraction_point,reflection_point,Tx];
                Incident_angle = angleOfIncidence(multipathtmp);
                ArrayOfPlanes = [1,reflecting_plane_eq];
                [reflectionLoss_dB,diffuseLoss,scatterPaths,tmp_plane,plane_mat] = getDbarcReflectionLoss(materialLibrary,... 
		                CADop(i, 14),multipathtmp,Incident_angle,mat_confg,IndoorSwitch, ...
                        polarization,Jones,scattering,ArrayOfPlanes,tmp_plane,plane_mat); 
                % Get diffraction loss (already calculated)
                %diffraction_loss_dB = path_data(1); % Assuming loss is the first column
        
                % Total path loss
                totalLoss = reflectionLoss_dB + diffraction_loss_dB;
                dod = Tx - reflection_point;
                doa = diffraction_point - Rx;
        
                % e. Construct the final output vector
                % This should match the 22-column format.
                % This is an example; adjust indices as per your final format.
                reconstructed_path = zeros(1, 22);
                reconstructed_path(1) = -2;
                % dod - direction of departure
                reconstructed_path(2:4) = dod;
                % doa - direction of arrival
                reconstructed_path(5:7) = doa;
                % Time delay
                reconstructed_path(8) = delay;
                
                %reconstructed_path(9) = 20*log10(wavelength / (4*pi*distance)) ...
                            %- reflectionLoss;            
                % Aod azimuth
                reconstructed_path(10)=atan2(dod(2),dod(1))*180/pi;
                % Aod elevation
                reconstructed_path(11)=atan2(dod(3),sqrt(dod(1)^2+dod(2)^2))*180/pi;
                % Aoa azimuth
                reconstructed_path(12)=atan2(doa(2),doa(1))*180/pi;
                % Aoa elevation
                reconstructed_path(13)=atan2(doa(3),sqrt(doa(1)^2+doa(2)^2))*180/pi;
    
                % Friis transmission loss
                 if ~isempty(LOS_output)% still assume LOS path exists
                     rel_aaod=abs(reconstructed_path(10)-LOS_AAOD);
                     rel_eaod=abs(reconstructed_path(11)-LOS_EAOD);
                     rel_aaoa=abs(reconstructed_path(12)-LOS_AAOA);
                     rel_eaoa=abs(reconstructed_path(13)-LOS_EAOA);
                     if strcmp(antenna, 'dir')
                         Gt=horn_gain(rel_eaod,rel_aaod);
                         Gr=horn_gain(rel_eaoa,rel_aaoa);
                     elseif strcmp(antenna, 'omni')
                         Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                         Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
                     elseif strcmp(antenna, 'scan')
                         Gt=23.7;
                         Gr=23.7;
                     end
                 elseif isempty(LOS_output) && (strcmp(antenna, 'dir') || strcmp(antenna, 'scan')) 
                     Gt=23.7;
                     Gr=23.7;
     % %If LOS path does not exist in directional channel 
     % sounding schemes, the Tx/Rx antennas will point on the 
     % best NLOS direction, but it's done in post processing
                     % LOS_AAOD=reconstructed_path(10);
                     % LOS_EAOD=reconstructed_path(11);
                     % LOS_AAOA=reconstructed_path(12);
                     % LOS_EAOA=reconstructed_path(13);
                     % LOS_output=[LOS_AAOD,LOS_EAOD,LOS_AAOA,LOS_EAOA];
                 elseif isempty(LOS_output) && strcmp(antenna, 'omni')
                     rel_aaod=abs(reconstructed_path(10)-pseudo_LOS_AAOD);
                     rel_eaod=abs(reconstructed_path(11)-pseudo_LOS_EAOD);
                     rel_aaoa=abs(reconstructed_path(12)-pseudo_LOS_AAOA);
                     rel_eaoa=abs(reconstructed_path(13)-pseudo_LOS_EAOA);
                     Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                     Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
                 end
    
                if strcmp(polarization, 'V-V') || strcmp(polarization, 'H-H')
                    reconstructed_path(9)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*total_distance)) - totalLoss;
                    reconstructed_path(16)=reconstructed_path(9) - Ptx - Gt - Gr;
                    Jones_dB = 10*log10(Jones);
                    reconstructed_path(17) = NaN;%cross-polarized channel gain,reflection
                    reconstructed_path(19) = min(Jones_dB + reconstructed_path(9));
                    reconstructed_path(22) = NaN;
                %elseif strcmp(polarization, 'cross') || strcmp(polarization, 'dual')
                elseif strcmp(polarization, 'dual')
                    %XPD = 34;
                    Jones_dB = 10*log10(Jones);
                    reconstructed_path(9)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*distance))-totalLoss(1);%co-polarized channel gain
                    reconstructed_path(16)=reconstructed_path(9)-Ptx-Gt-Gr;
                    reconstructed_path(17)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*distance))-totalLoss(2);%another co-polarization channel gain
                    reconstructed_path(19)= min(Jones_dB(:,1) + reconstructed_path(9));% cross-polarization channel gain at Tx
                    reconstructed_path(22)= min(Jones_dB(:,1) + reconstructed_path(17)); % cross-polarization channel gain at Rx
                end
                % if reconstructed_path(17)<=-120
                %     reconstructed_path(17)=NaN;
                % end
                if reconstructed_path(19)<=-120
                    reconstructed_path(19)=NaN;
                end
                if reconstructed_path(22)<=-120
                    reconstructed_path(22)=NaN;
                end
                if enablePhase
                    reconstructed_path(18) = mod(1*pi+total_distance/wavelength*2*pi,2*pi);
                else
                    reconstructed_path(18) = 0;
                end
                reconstructed_path(14)=Gt;
                reconstructed_path(15)=Gr;
                reconstructed_path(20) = dopplerFactor * frequency;
                reconstructed_path(21) = 0;
                if reconstructed_path(9)<=-120 && (reconstructed_path(17)<=-120 || isnan(reconstructed_path(17)))
                    reconstructed_path=[];
                    multipath1 = [];
                    % output(indexOutput).reconstructed_path=reconstructed_path;
                    % multipath(iterateNumberOfRowsArraysOfPlanes,:)=[];
                else
                    if reconstructed_path(17)<=-120
                        reconstructed_path(17)=NaN;
                    end
                    if reconstructed_path(9)<=-120
                        reconstructed_path(9)=NaN;
                    end
                    output = [output; reconstructed_path];
                    multipath1 = [-2,Rx,diffraction_point,reflection_point,Tx];
                    multipath = [multipath;multipath1];
                end
                % If reconstructed_path detectable, generate compone
                % You can also populate 'multipath' if needed
            end
        end
    end
end

function incidentAngle = angleOfIncidence(multipath) 
% ANGLEOFINCIDENCE returns angle of incident for first and second order
% reflections
% 
% Inputs:
% multipath - consists of specular multipath parameters. This vector is
% used to calculate angle(s) of incident.
% 
% Output: 
% incidentAngle - incident angles till second order reflections

differenceVectorRxFirstPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
% differenceVectorRxFirstPoI is the difference vector between Rx 
% and first point of incidence (PoI).
differenceVectorTxFirstPoI = (multipath(1,8:10) - multipath(1,5:7))...
                            /( norm(multipath(1,8:10) - multipath(1,5:7)));
% differenceVectorTxFirstPoI is the difference vector between Tx 
% and first point of incidence (PoI). 
dpAoI = dot(differenceVectorRxFirstPoI, differenceVectorTxFirstPoI);
% dpAoI is the dot product between differenceVectorRxFirstPoI and 
% differenceVectorTxFirstPoI will give the cos of angle between the vectors.
incidentAngle(1) = 0.5*acosd(dpAoI);
% Half of the angle between the vectors differenceVectorRxFirstPoI and 
% differenceVectorTxFirstPoI is the angle of incidence. This is because 
% angle of incidence is equal to angle of reflection.
if multipath(1,1) == 2 % This is for second order reflection
    differenceVectorRxSecondPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
    differenceVectorFirstPoISecondPoI =...  
                            (multipath(1,8:10) - multipath(1,5:7))...
                            /( norm(multipath(1,8:10) - multipath(1,5:7)));
    differenceVectorSecondPoIFirstPoI = -differenceVectorFirstPoISecondPoI;
    differenceVectorTxFirstPoI =...
        (multipath(1,11:13) - multipath(1,8:10))/...
        norm((multipath(1,11:13) - multipath(1,8:10)));
    dpAoI = dot(differenceVectorSecondPoIFirstPoI,...
                differenceVectorTxFirstPoI);
    incidentAngle(1) = 0.5*acosd(dpAoI);
    dpAoI = dot(differenceVectorRxSecondPoI,...
                differenceVectorFirstPoISecondPoI);
    incidentAngle(2)  = 0.5*acosd(dpAoI);
end
if multipath(1,1) == 3 % This is for third order reflection
    differenceVectorRxThirdPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
    differenceVectorRxSecondPoI = (multipath(1,5:7) - multipath(1,8:10))...
                            /( norm(multipath(1,5:7) - multipath(1,8:10)));
    differenceVectorFirstPoISecondPoI = (multipath(1,8:10) - multipath(1,11:13))...
                            /( norm(multipath(1,8:10) - multipath(1,11:13)));
    differenceVectorTxFirstPoI =...
        -(multipath(1,11:13) - multipath(1,14:16))/...
        norm((multipath(1,11:13) - multipath(1,14:16)));
    dpAoI = dot(differenceVectorTxFirstPoI,...
                differenceVectorFirstPoISecondPoI);
    incidentAngle(1) = 0.5*acosd(dpAoI);
    dpAoI = dot(-differenceVectorFirstPoISecondPoI,...
                differenceVectorRxSecondPoI);
    incidentAngle(2)  = 0.5*acosd(dpAoI);
    dpAoI = dot(-differenceVectorRxSecondPoI,...
                differenceVectorRxThirdPoI);
    incidentAngle(3)  = 0.5*acosd(dpAoI);

end
end