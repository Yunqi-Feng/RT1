function [CADOutput, materialSwitch,scatterer_confg] = xmlreader(filename, ...
    MaterialLibrary, referencePoint, r, IndoorSwitch,confg_dim,...
    scatterer_num,scatterer_dim,scatterer_dis,NodeLocTx, ...
    NodeLocRx,scatterer_confg,scatterer_height_max,mat_confg,...
    antenna,Ptx,polarization,scattering)
% XMLREADER function extracts the information of CAD file (AMF). 
% 
% Inputs:
% filename - file name 
% MaterialLibrary - material database with all the material parameters
% referencePoint - center of the sphere 
% r - radius of the sphere
% IndoorSwitch - defines whether the scenario is indoor or not
% confg_dim:if indoor, room dimension. if outdoor, reference dimension
% s: the floor map which has been converted to struct format
% x_cen/y_cen/z_cen: the centroid of scatterers
% x_dep/y_dep/z_dep: the size of scatterers
% scatterer_num: scatterer number
% scatterer_dim: scatterer dimension, 0=small, 1=medium, 2=big
% scatterer_dis: scatterer distribution
% NodeLocTx: AP location(s)
% NodeLocRx: UE location(s)

% Outputs:
% CADOutput - contains all the extracted triangles
% materialSwitch - a boolean to know whether the material information is
% present in the CAD file. If any one of the materials is missing 
% XMLREADER function returns materialSwitch = 0
% scatterer_confg: configured scatterers in the floor map described by
% their centers and dimensions


%--------------------------Software Disclaimer-----------------------------
%
% NIST-developed software is provided by NIST as a public service. You may 
% use, copy and distribute copies of the software in any medium, provided 
% that you keep intact this entire notice. You may improve, modify and  
% create derivative works of the software or any portion of the software, 
% and you  may copy and distribute such modifications or works. Modified 
% works should carry a notice stating that you changed the software and  
% should note the date and nature of any such change. Please explicitly  
% acknowledge the National Institute of Standards and Technology as the 
% source of the software.
% 
% NIST-developed software is expressly provided "AS IS." NIST MAKES NO
% WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION  
% OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND 
% DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF 
% THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS 
% WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS  
% REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT 
% NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF 
% THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with  
% its use, including but not limited to the risks and costs of program 
% errors, compliance with applicable laws, damage to or loss of data, 
% programs or equipment, and the unavailability or interruption of 
% operation. This software is not intended to be used in any situation  
% where a failure could cause risk of injury or damage to property. The 
% software developed by NIST employees is not subject to copyright 
% protection within the United States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Used MATLAB functions 
%   instead of custom ones, improved MaterialLibrary access, readibility, 
%   performance in general
% Modified by: Neeraj Varshney <neeraj.varshney@nist.gov>, support multiple
%   objects and different length units in amf file
% Modified by: Yunqi Feng <Yunqi.Feng@UGent.be>, support flexible
% configuration
% December 2023
s = xml2struct(filename);
s.amf.metadata.Text='Yunqi Feng';
if scatterer_num > 0
    if isempty(scatterer_confg)
        % Call the new robust placement function
        [x_cen, y_cen, z_cen, x_dep, y_dep, z_dep] = generate_valid_scatterer_layout(scatterer_num, confg_dim, scatterer_dim, scatterer_dis, NodeLocTx, NodeLocRx, IndoorSwitch, scatterer_confg, scatterer_height_max);

        % Update the final scatterer number (should be the same unless an error was thrown)
        scatterer_num = length(x_cen);
        
        % Generate the geometry for the validated scatterers
        s = scatterer_generation(s, confg_dim, scatterer_num, x_cen, y_cen, z_cen, x_dep, y_dep, z_dep, IndoorSwitch);
        
        % Store the final valid configuration
        scatterer_confg = [x_cen; y_cen; z_cen; x_dep; y_dep; z_dep];
    else
        % Logic to handle pre-configured scatterers remains the same
        x_cen=scatterer_confg(1,1:scatterer_num);
        y_cen=scatterer_confg(2,1:scatterer_num);
        z_cen=scatterer_confg(3,1:scatterer_num);
        x_dep=scatterer_confg(4,1:scatterer_num);
        y_dep=scatterer_confg(5,1:scatterer_num);
        z_dep=scatterer_confg(6,1:scatterer_num);
    end
else
    s = scenario_generation(s, confg_dim, IndoorSwitch);
    scatterer_confg = [];
end
% if scatterer_num>0
%     [x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,scatterer_num]=scatterer_config(confg_dim,scatterer_num,...
%     scatterer_dim,scatterer_dis,NodeLocTx,NodeLocRx,IndoorSwitch, ...
%     scatterer_confg,scatterer_height_max);
%     scatterer_num=length(x_cen);
%     s=scatterer_generation(s,confg_dim,scatterer_num,...
%     x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,IndoorSwitch);
%     % for i=1:scatterer_num
%     %     idx=1;
%     %     while (NodeLocTx(1)>x_cen(i)-0.5*x_dep(i)&&NodeLocTx(1)<x_cen(i)+0.5*x_dep(i)&&...
%     %             NodeLocTx(2)>y_cen(i)-0.5*y_dep(i)&&NodeLocTx(2)<y_cen(i)+0.5*y_dep(i)&&...
%     %             NodeLocTx(3)>z_cen(i)-0.5*z_dep(i)&&NodeLocTx(3)<z_cen(i)+0.5*z_dep(i))||...
%     %             (NodeLocRx(1)>x_cen(i)-0.5*x_dep(i)&&NodeLocRx(1)<x_cen(i)+0.5*x_dep(i)&&...
%     %             NodeLocRx(2)>y_cen(i)-0.5*y_dep(i)&&NodeLocRx(2)<y_cen(i)+0.5*y_dep(i)&&...
%     %             NodeLocRx(3)>z_cen(i)-0.5*z_dep(i)&&NodeLocRx(3)<z_cen(i)+0.5*z_dep(i))
%     %         [x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,scatterer_num]=scatterer_config(confg_dim,scatterer_num,...
%     %     scatterer_dim,scatterer_dis,NodeLocTx,NodeLocRx,IndoorSwitch, ...
%     %     scatterer_confg,scatterer_height_max);
%     %         s=scatterer_generation(s,confg_dim,scatterer_num,...
%     %     x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,IndoorSwitch);
%     %         idx=idx+1;
%     %         if idx>10
%     %             error('inappropriate configuration, reduce the number or dimension of scatterers')
%     %         end
%     %     end 
%     % end
%     % In xmlreader.m, a potential modification
%     % for i=1:scatterer_num
%     %     is_overlapping = true;
%     %     attempt = 0;
%     % 
%     %     while is_overlapping && attempt < 20 % Try 20 times for just this one scatterer
%     % 
%     %         % Bounding box check for the i-th scatterer
%     %         tx_in = (NodeLocTx(1) > x_cen(i)-0.5*x_dep(i) && NodeLocTx(2) > y_cen(i)-0.5*y_dep(i) && NodeLocTx(3) > z_cen(i)-0.5*z_dep(i)); % (shortened for clarity)
%     %         rx_in = (NodeLocRx(1) > x_cen(i)-0.5*x_dep(i) && NodeLocRx(2) > y_cen(i)-0.5*y_dep(i) && NodeLocRx(3) > z_cen(i)-0.5*z_dep(i)); % (shortened for clarity)
%     % 
%     %         if tx_in || rx_in
%     %             % Regenerate ONLY the properties for the i-th scatterer
%     %             [new_xc, new_yc, new_zc, new_xd, new_yd, new_zd, ~] = scatterer_config(confg_dim,scatterer_num,scatterer_dim,scatterer_dis, ...
%     %                 NodeLocTx,NodeLocRx,IndoorSwitch,scatterer_confg,scatterer_height_max);
%     % 
%     %             % Update the properties for just this scatterer
%     %             % x_cen(i) = new_xc; y_cen(i) = new_yc; z_cen(i) = new_zc;
%     %             % x_dep(i) = new_xd; y_dep(i) = new_yd; z_dep(i) = new_zd;
%     %             x_cen(i) = new_xc(i); y_cen(i) = new_yc(i); z_cen(i) = new_zc(i);
%     %             x_dep(i) = new_xd(i); y_dep(i) = new_yd(i); z_dep(i) = new_zd(i);
%     % 
%     %             attempt = attempt + 1;
%     %             is_overlapping = true;
%     %         else
%     %             is_overlapping = false;
%     %         end
%     %     end
%     % 
%     %     if attempt >= 100
%     %         error('Could not place scatterer %d without overlapping Tx/Rx.', i);
%     %     end
%     indices_to_delete = [];
% 
%     % Loop through each finalized scatterer to check for overlap
%     for i = 1:scatterer_num
%         % Bounding Box Check for scatterer 'i'
%         tx_is_inside = (NodeLocTx(1) >= x_cen(i) - 0.5 * x_dep(i) && NodeLocTx(1) <= x_cen(i) + 0.5 * x_dep(i) && ...
%                         NodeLocTx(2) >= y_cen(i) - 0.5 * y_dep(i) && NodeLocTx(2) <= y_cen(i) + 0.5 * y_dep(i) && ...
%                         NodeLocTx(3) >= z_cen(i) - 0.5 * z_dep(i) && NodeLocTx(3) <= z_cen(i) + 0.5 * z_dep(i));
% 
%         rx_is_inside = (NodeLocRx(1) >= x_cen(i) - 0.5 * x_dep(i) && NodeLocRx(1) <= x_cen(i) + 0.5 * x_dep(i) && ...
%                         NodeLocRx(2) >= y_cen(i) - 0.5 * y_dep(i) && NodeLocRx(2) <= y_cen(i) + 0.5 * y_dep(i) && ...
%                         NodeLocRx(3) >= z_cen(i) - 0.5 * z_dep(i) && NodeLocRx(3) <= z_cen(i) + 0.5 * z_dep(i));
% 
%         if tx_is_inside || rx_is_inside
%             fprintf('    WARNING: Scatterer %d encloses a Tx/Rx node. Marking for removal.\n', i);
%             indices_to_delete = [indices_to_delete, i];
%         end
%     end
% 
%     % If any scatterers were marked for deletion, remove them
%     if ~isempty(indices_to_delete)
%         % Remove the conflicting scatterers from the properties arrays
%         x_cen(indices_to_delete) = [];
%         y_cen(indices_to_delete) = [];
%         z_cen(indices_to_delete) = [];
%         x_dep(indices_to_delete) = [];
%         y_dep(indices_to_delete) = [];
%         z_dep(indices_to_delete) = [];
% 
%         % Update the official scatterer count
%         num_deleted = length(indices_to_delete);
%         scatterer_num = scatterer_num - num_deleted;
%         warning('Removed %d overlapping scatterer(s). Final scatterer count is now %d.', num_deleted, scatterer_num);
% 
%         % The cleanest way to update the scene is to reload the base environment
%         % and regenerate ONLY the valid scatterers.
%         % s = xml2struct(filename); % Reload the original file
%         % s.amf.metadata.Text='Yunqi Feng';
%         % if isempty(scatterer_confg)
%         %      s = scenario_generation(s,confg_dim,IndoorSwitch); % Re-add walls/floor
%         %      % s = scatterer_generation(s,confg_dim,scatterer_num, ...
%         %      %     x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,IndoorSwitch);
%         % end
%     end
% else
%     s=scenario_generation(s,confg_dim,IndoorSwitch);
%     scatterer_confg=[];
% end
        
        % After the loop, all scatterers have valid positions.
        % Now, generate the geometry for all of them at once.
%         s = scatterer_generation(s, confg_dim, scatterer_num, x_cen, y_cen, z_cen, x_dep, y_dep, z_dep, IndoorSwitch);
%         scatterer_confg_tmp=[x_cen;y_cen;z_cen;x_dep;y_dep;z_dep];
%         if ~isempty(scatterer_confg)
%             scatterer_confg=scatterer_confg;
%         else
%             scatterer_confg=[scatterer_confg scatterer_confg_tmp];
%         end
%     end
% else
%     s=scenario_generation(s,confg_dim,IndoorSwitch);
%     scatterer_confg=[];
% end
if isfield(s.amf, 'material')
    materialSwitch = 1;
    
    sizeMaterials1 = size(s.amf.material');
    if sizeMaterials1(2)>1 &&  sizeMaterials1(1)==1
        sizeMaterials = sizeMaterials1;
    else
        sizeMaterials = sizeMaterials1(1);
    end
    
else
    materialSwitch = 0;
    
end
%% Iterating through all the subdivisions (volumes) and extracting the triangle information

CADOutput = [];
lengthObject = length(s.amf.object);
if lengthObject~=0
    for iterateObjects = 1:lengthObject                            % For multiple objects
        if lengthObject == 1
            volume = s.amf.object{1,1}.mesh.volume';
            sizeVolume = size(volume);
            sObject = s.amf.object{1,1};
        else
            volume = s.amf.object{1, iterateObjects}.mesh.volume';
            sizeVolume = size(volume);
            sObject = s.amf.object{1, iterateObjects};
        end
        for iterateVolume = 1:sizeVolume
            if sizeVolume(1) ~= 1
                triangles = sObject.mesh.volume{1, iterateVolume}.triangle';
            else
                triangles = sObject.mesh.volume.triangle';
            end
        
            if materialSwitch == 1
                if sizeVolume(1) ~= 1
                    materialId = sObject.mesh.volume{1, iterateVolume}.Attributes.materialid;
                else
                    materialId = sObject.mesh.volume.Attributes.materialid;
                end
            
                for iterateMaterials = 1:sizeMaterials
                    if sizeMaterials ~= 1
                        if str2double(materialId) == str2double...
                                (s.amf.material{1, iterateMaterials}.Attributes.id)
                            material = s.amf.material{1, iterateMaterials}.metadata.Text;
                        end
                    
                    elseif sizeVolume(1) == 1 && sizeMaterials == 1

                        if str2double(materialId) == str2double(s.amf.material{1,1}.Attributes.id)
                            material = s.amf.material{1,1}.metadata.Text;
                        end
                    
                    end
                end
            end
        %% Extracting the vertices information of the triangles
        
            sizeTriangle = size(triangles);
            CADOutputTemp = [];
            for iterateTriangles = 1:sizeTriangle(1)
                if isfield(s.amf, 'Attributes')
                    switch s.amf.Attributes.unit
                        case 'micrometer'
                            unitConversion = 1e-6;
                        case 'millimeter'
                            unitConversion = 1e-3;
                        case 'meter'
                            unitConversion = 1;
                        case 'kilometer'
                            unitConversion = 1e3;
                        case 'inch'
                            unitConversion = 0.0254;
                        case 'foot'
                            unitConversion = 0.3048;
                        case 'mile'
                            unitConversion = 1609.34;
                        otherwise
                            error('xmlreader does not support this unit.');
                    end
                    v1 = getTriangleVertex(sObject, iterateVolume, ...
                    iterateTriangles, 'v1', sizeVolume)*unitConversion;
                    v2 = getTriangleVertex(sObject, iterateVolume, ...
                    iterateTriangles, 'v2', sizeVolume)*unitConversion;
                    v3 = getTriangleVertex(sObject, iterateVolume, ...
                    iterateTriangles, 'v3', sizeVolume)*unitConversion;
                else
                    error('Length unit is missing in the xml/amf file. Add <amf  unit="?" in Line 1>');
                end
                
            
            % Calculating the plane equation of triangles
            
                vector1 = v2 - v3;
                vector2 = -(v2 - v1);
            
                normal = cross(vector2, vector1) * (1-(2*IndoorSwitch));
                normal = round(normal/norm(normal), 4);
                vector3 = v2;
            % for box. remove for others
                D = -dot(normal, vector3);
            
            % Storing Material information in output if the material exists in the material database
                if materialSwitch==1
                    materialFound = false;
                %material='Scatterer';
                
                    for iterateMaterials=1:size(MaterialLibrary, 1)
                        if strcmpi(MaterialLibrary.Reflector{iterateMaterials}, material)
                            CADOutputTemp(14) = iterateMaterials;
                            materialFound = true;
                            break
                        end
                    end
                
                % Storing triangle vertices and plane equations in output
                % Part where output file is created. It contains the triangle vertices
                % in first nine columns, plane equations in the next four columns
                    if ~materialFound
                        materialSwitch=0;
                        warning('Material ''%s'' not found. Disabling materials', material)
                    end
                
                end
            
            %
                if materialSwitch == 0 && size(CADOutput, 2) == 14
                    CADOutput(:, 14) = [];
                end
            
                CADOutputTemp(1:3) = round(v1, 6);
                CADOutputTemp(4:6) = round(v2, 6);
                CADOutputTemp(7:9) = round(v3, 6);
                CADOutputTemp(10:12) = round(normal, 4);
                CADOutputTemp(13) = round(D, 4);
            
            % We are using distance limitation at this step
                if isinf(r)
                    [switchDistance] = 1;
                else
                    [switchDistance] = verifydistance(r, referencePoint, CADOutputTemp, 1);
                end
            
            % If the triangles are within the given distance we increase the count,
            % else the next triangle will replace the present row (as count remains constant)
                if switchDistance==1
                    CADOutput = [CADOutput; CADOutputTemp];
                end
            
            end
        end
    end
end
end


%% Utils
function v = getTriangleVertex(sObject, volumeIdx, triangIdx, vertexIdx, sizeVolume)

if sizeVolume(1) ~= 1
    vertex = str2double(sObject.mesh.volume{1, volumeIdx}.triangle{1, triangIdx}.(vertexIdx).Text)+1;
else
    vertex = str2double(sObject.mesh.volume.triangle{1, triangIdx}.(vertexIdx).Text)+1;
end

x = str2double(sObject.mesh.vertices.vertex{1, vertex}.coordinates.x.Text);
y = str2double(sObject.mesh.vertices.vertex{1, vertex}.coordinates.y.Text);
z = str2double(sObject.mesh.vertices.vertex{1, vertex}.coordinates.z.Text);
v = [x, y, z];

end