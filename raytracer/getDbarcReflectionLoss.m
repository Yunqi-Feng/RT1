function [reflectionLoss,diffuseLoss,scatterPaths,tmp_plane,plane_mat] = getDbarcReflectionLoss(materialLibrary, arrayOfMaterials, ...
    multipath,Incident_angle,mat_confg,IndoorSwitch,polarization,...
    Jones,scattering,ArrayOfPlanes,tmp_plane,plane_mat)
% GETDBARCREFLECTIONLOSS returns the ray reflection loss based on IMEC 
% D-band radio measurements
% 
% reflectionLoss = GETDBARCREFLECTIONLOSS(materialLibrary, arrayOfMaterials,multipath)
% These are my matlab codes to compute the reflection loss, which is up to second order
% reflection. Some useful parameters are: 1. materialLibrary, which contains the refractive
% index, absorption rate and roughness of materials. The parameter for the walls,
% the floors and the ceiling are set to 0 which can be convenient for the following assignment. 
% 2. arrayOfMaterials, which are the material index from the library of interactive  surfaces. 
% Currently second order reflection can be modeled, so its length is up to two. 
% 3. Incident Anlge, which is one angle for the first order reflection and two angles 
% for second order reflection. 4. mat_confg: refer to names to the material and will be 
% mapped to exact material parameters by the switch command in line 32. 5. polarization 
% and Jones: determine the reflection loss of different wave polarizaztions. 
% 6. arrayofPlanes, the first element is the reflection order, the subsequent four elements 
% are the interaction plane,e.g., 0 0 1 1 means 0*x+0*y+z=1, if first order reflection is 
% calculated. There are nine elements in the arrayofPlanes if second order reflections are 
% calculated, following the above description. 7. tmp_planes: although scatterers in the 
% environment are all named 'scatterers', different scatterers should be assigned with 
% different materials, and the previous scatterer involved in the interaction should maintain 
% the same material properties. So tmp_plame/plane mat is used to check this problem, which stores 
% the coordinate and materials of the scatterer plane respectively. 8 scattering are not considered 
% in the modeling so just negelect them. Now I want to furthur model the third order reflection. 
% Following my description, extend the matlab codes to model third order reflections.






    %% Init
u0=4*pi*10^(-7);
c=3e8;
f=140e9;
lab=c/f;
e0=1/(u0*c^2);
Z0=377;
reflectionLoss = 0;
orderReflection=multipath(1,1);
reflectionCoefficient=ones(2,orderReflection);
% if orderReflection==1
%     normal = ArrayOfPlanes(2:5);
% elseif orderReflection==2
%     normal1 = ArrayOfPlanes(2:5);
%     normal2 = ArrayOfPlanes(end-3:end);
%     normal = [normal1;normal2];
% end
normal = ArrayOfPlanes(end-3:end-1);
if IndoorSwitch==1
    mat_indoor=[1,2,3,4,5,6];%index in the excel for material assignment
    for i=1:length(mat_indoor)
        if materialLibrary.permittivity(mat_indoor(i))==0
            switch mat_confg(i)% different kinds of materials
                case 1
                    tmp=17;
                case 2
                    tmp=randi([15 16]);
                case 3
                    tmp=18;
                case 4
                    index=[8 9 10 11 12 13 14];
                    tmp=index(randperm(numel(index),1));
                otherwise
                    error('Undefined material')
            end
            materialLibrary.permittivity(mat_indoor(i))=materialLibrary.permittivity(tmp);
            materialLibrary.conductivity(mat_indoor(i))=materialLibrary.conductivity(tmp);
            materialLibrary.roughness(mat_indoor(i))=materialLibrary.roughness(tmp);
        end
    end
    if length(arrayOfMaterials)==1
        if arrayOfMaterials==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        end
    elseif length(arrayOfMaterials)==2
        if arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(2)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(2)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1  && arrayOfMaterials(1)~=7
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1)
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,scatterer2_index]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(scatterer2_index,1);
            scatterer2_con=plane_mat(scatterer2_index,2);
            scatterer2_rou=plane_mat(scatterer2_index,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1)
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1)
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows'); 
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        end
    elseif length(arrayOfMaterials)==3
        if arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(2)~=7 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(2)~=7 && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1  && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,scatterer2_index]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(scatterer2_index,1);
            scatterer2_con=plane_mat(scatterer2_index,2);
            scatterer2_rou=plane_mat(scatterer2_index,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows'); 
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)~=7 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)~=7 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        end
    end
else
    mat_indoor=6;
    for i=1:length(mat_indoor)
        if materialLibrary.permittivity(mat_indoor(i))==0
            materialLibrary.permittivity(mat_indoor(i))=3;
            materialLibrary.conductivity(mat_indoor(i))=1;
            materialLibrary.roughness(mat_indoor(i))=0.1*rand;
        end
    end
    if length(arrayOfMaterials)==1
        if arrayOfMaterials==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        end
    elseif length(arrayOfMaterials)==2
        if arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(2)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(2)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1  && arrayOfMaterials(1)~=7
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1)
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,scatterer2_index]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(scatterer2_index,1);
            scatterer2_con=plane_mat(scatterer2_index,2);
            scatterer2_rou=plane_mat(scatterer2_index,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1)
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1)
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows'); 
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        end
    elseif length(arrayOfMaterials)==3
        if arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(2)~=7 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(2)~=7 && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1  && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,scatterer2_index]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(scatterer2_index,1);
            scatterer2_con=plane_mat(scatterer2_index,2);
            scatterer2_rou=plane_mat(scatterer2_index,3);
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(3)~=7
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif arrayOfMaterials(1)==7 && arrayOfMaterials(2)==7 && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1) && (~isempty(tmp_plane) && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1) && arrayOfMaterials(3)~=7
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows'); 
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer2_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer2_index);
            scatterer2_con=materialLibrary.conductivity(scatterer2_index);
            scatterer2_rou=materialLibrary.roughness(scatterer2_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)~=7 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(2)~=7 && arrayOfMaterials(1)~=7 && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane,'rows'))==1 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            [~,idx]= ismember(ArrayOfPlanes(2:5),tmp_plane, 'rows');
            materialLibrary.permittivity(7)=plane_mat(idx,1);
            materialLibrary.conductivity(7)=plane_mat(idx,2);
            materialLibrary.roughness(7)=plane_mat(idx,3);
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            [~,idx]= ismember(ArrayOfPlanes(6:9),tmp_plane, 'rows');
            scatterer2_per=plane_mat(idx,1);
            scatterer2_con=plane_mat(idx,2);
            scatterer2_rou=plane_mat(idx,3);
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane,'rows'))==1 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            [~,idx]= ismember(ArrayOfPlanes(10:13),tmp_plane, 'rows');
            scatterer3_per=plane_mat(idx,1);
            scatterer3_con=plane_mat(idx,2);
            scatterer3_rou=plane_mat(idx,3);
        elseif ~isempty(tmp_plane) && arrayOfMaterials(3)==7 && sum(ismember(ArrayOfPlanes(10:13),tmp_plane))~=4 && arrayOfMaterials(2)==7 && sum(ismember(ArrayOfPlanes(2:5),tmp_plane))~=4 && arrayOfMaterials(1)==7 && sum(ismember(ArrayOfPlanes(6:9),tmp_plane))~=4
            index=[10,12,13,14];
            scatterer_index=index(randperm(numel(index),1));
            materialLibrary.permittivity(7)=materialLibrary.permittivity(scatterer_index);
            materialLibrary.conductivity(7)=materialLibrary.conductivity(scatterer_index);
            materialLibrary.roughness(7)=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(2:5)];
            plane_mat=[plane_mat;materialLibrary.permittivity(7) materialLibrary.conductivity(7) materialLibrary.roughness(7)];
            scatterer_index=index(randperm(numel(index),1));
            scatterer2_per=materialLibrary.permittivity(scatterer_index);
            scatterer2_con=materialLibrary.conductivity(scatterer_index);
            scatterer2_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(6:9)];
            plane_mat=[plane_mat;scatterer2_per scatterer2_con scatterer2_rou];
            scatterer_index=index(randperm(numel(index),1));
            scatterer3_per=materialLibrary.permittivity(scatterer_index);
            scatterer3_con=materialLibrary.conductivity(scatterer_index);
            scatterer3_rou=materialLibrary.roughness(scatterer_index);
            tmp_plane=[tmp_plane;ArrayOfPlanes(10:13)];
            plane_mat=[plane_mat;scatterer3_per scatterer3_con scatterer3_rou];
        end
    end
end
Rou=[];
%% Loop over reflection order
for reflectionOrderIndex = 1:orderReflection
    if reflectionOrderIndex == 1
        reflectionMaterialIndex = arrayOfMaterials(1,1);
    elseif  reflectionOrderIndex == 2
        reflectionMaterialIndex = arrayOfMaterials(1,2);
    elseif reflectionOrderIndex == 3
        reflectionMaterialIndex = arrayOfMaterials(1,3);
    end
    if (reflectionMaterialIndex==7 && (exist('scatterer2_per','var') || exist('scatterer3_per','var')) && (reflectionOrderIndex==2 || reflectionOrderIndex==3))
        if ~(exist('scatterer3_per','var') && reflectionOrderIndex==3)
            per = scatterer2_per;
            con = scatterer2_con;
            per_img = con / (2 * pi * f * e0);
            n = sqrt((sqrt(per^2+per_img^2) + per)/2);
            k = sqrt((sqrt(per^2+per_img^2) - per)/2);
            a = 4 * pi *f / c * k;
            % n=scatterer2_per;
            % absorption_coeff=scatterer2_con;
            % a=absorption_coeff*100;
            % roughness=scatterer2_rou;
            % r=roughness*10^(-3);
            aor = Incident_angle(reflectionOrderIndex);
            %Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(2*pi*f))));
            Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(4*pi*f))));
            aot = asind(sind(aor)*Z/Z0);
            %rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2)*besselj(0,8*(pi*r*cos(aor)/lab)); % Beckmann-Kirchhoff theory in the reflection direction      
            rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2);
            Rou = [Rou,rou];
            reflectionCoefficient(:, reflectionOrderIndex) = rou.*[ ...
            (Z*cosd(aor)-Z0*cosd(aot))/(Z*cosd(aor)+Z0*cosd(aot)); ... % Parallel 
            (Z*cosd(aot)-Z0*cosd(aor))/(Z*cosd(aot)+Z0*cosd(aor))];   % Perpendicular
        else
            per = scatterer3_per;
            con = scatterer3_con;
            per_img = con / (2 * pi * f * e0);
            n = sqrt((sqrt(per^2+per_img^2) + per)/2);
            k = sqrt((sqrt(per^2+per_img^2) - per)/2);
            a = 4 * pi *f / c * k;
            % n=scatterer3_per;
            % absorption_coeff=scatterer3_con;
            % a=absorption_coeff*100;
            % roughness=scatterer3_rou;
            r=roughness*10^(-3);
            aor = Incident_angle(reflectionOrderIndex);
            %Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(2*pi*f))));
            Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(4*pi*f))));
            aot = asind(sind(aor)*Z/Z0);
            %rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2)*besselj(0,8*(pi*r*cos(aor)/lab)); % Beckmann-Kirchhoff theory in the reflection direction      
            rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2);
            Rou = [Rou,rou];
            reflectionCoefficient(:, reflectionOrderIndex) = rou.*[ ...
            (Z*cosd(aor)-Z0*cosd(aot))/(Z*cosd(aor)+Z0*cosd(aot)); ... % Parallel 
            (Z*cosd(aot)-Z0*cosd(aor))/(Z*cosd(aot)+Z0*cosd(aor))];   % Perpendicular
        end
    else
        per = materialLibrary.permittivity(reflectionMaterialIndex);
        con = materialLibrary.conductivity(reflectionMaterialIndex);
        per_img = con / (2 * pi * f * e0);
        n = sqrt((sqrt(per^2+per_img^2) + per)/2);
        k = sqrt((sqrt(per^2+per_img^2) - per)/2);
        a = 4 * pi *f / c * k;
        % refractive_index = materialLibrary.permittivity(reflectionMaterialIndex); 
        % n=refractive_index;
        % absorption_coeff = materialLibrary.conductivity(reflectionMaterialIndex);
        % a=absorption_coeff*100;
        roughness=materialLibrary.roughness(reflectionMaterialIndex);
        r=roughness*10^(-3);
        aor = Incident_angle(reflectionOrderIndex);
        %Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(2*pi*f))));
        Z=sqrt(u0/(e0*(n^2-(a*c/(4*pi*f))^2-1j*2*n*a*c/(4*pi*f))));
        aot = asind(sind(aor)*Z/Z0);
        %rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2)*besselj(0,8*(pi*r*cos(aor)/lab));
        rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2);
        %rou=1;
        Rou = [Rou,rou];
        reflectionCoefficient(:, reflectionOrderIndex) = rou.*[ ...
        (Z*cosd(aor)-Z0*cosd(aot))/(Z*cosd(aor)+Z0*cosd(aot)); ... % Parallel 
        (Z*cosd(aot)-Z0*cosd(aor))/(Z*cosd(aot)+Z0*cosd(aor))];   % Perpendicular
    end
end
if orderReflection==1
    norm1 = ArrayOfPlanes(2:4);%normal vector of the interacted plane, dod the incident vector
    dod1 = (multipath(1,8:10) - multipath(1,5:7));
    mynorm = norm1;
    dod = dod1;
elseif orderReflection==2
    norm1 = ArrayOfPlanes(2:4);
    norm2 = ArrayOfPlanes(6:8);
    mynorm = [norm1;norm2];
    dod1 = (multipath(1,11:13) - multipath(1,8:10));
    dod2 = (multipath(1,8:10) - multipath(1,5:7));
    dod = [dod1;dod2];
elseif orderReflection==3
    norm1 = ArrayOfPlanes(2:4);
    norm2 = ArrayOfPlanes(6:8);
    norm3 = ArrayOfPlanes(10:12);
    mynorm = [norm1;norm2;norm3];
    dod1 = (multipath(1,14:16) - multipath(1,11:13));
    dod2 = (multipath(1,11:13) - multipath(1,8:10));
    dod3 = (multipath(1,8:10) - multipath(1,5:7));
    dod = [dod1;dod2;dod3];
end
inc_plane = [];
for i = 1:orderReflection
    inc_plane_tmp = cross(mynorm(i,:),dod(i,:));
    inc_plane = [inc_plane;inc_plane_tmp];
end
R = [];
if scattering == 0
    %1 - rou
    %rou=exp(-0.5*(4*pi*r*cosd(aor)/lab)^2)*besselj(0,8*(pi*r*cos(aor)/lab));% Beckmann-Kirchhoff theory in the reflection direction
    %scatteringcoefficient = 1 - rou;
    %reflectionCoefficient = prod(sqrt(abs(reflectionCoefficient).^2),2);
    if strcmp(polarization,'V-V')
        E_inc = [0,0,1];
        for i = 1:orderReflection
            is_perpendicular = norm(cross(E_inc, inc_plane(i,:)));
            if is_perpendicular < 1e-10
                R = [R,reflectionCoefficient(2,i)];
            else
                R = [R,reflectionCoefficient(1,i)];
            end
        end
        reflectionCoefficient = prod(sqrt(abs(R).^2));
        %R = diag(reflectionCoefficient) * Jones;
        %reflectionLoss = sum(R(:));
        %reflectionCoefficient = reflectionCoefficient(1,:);
        %reflectionCoefficient = diag(reflectionCoefficient) * Jones;
        %reflectionLoss = prod(abs(reflectionCoefficient).^2); 
        reflectionLoss = -20*log10(reflectionCoefficient);
    if strcmp(polarization,'H-H')
        % R = diag(reflectionCoefficient) * Jones;
        % reflectionLoss = sum(R(:));
        % %reflectionCoefficient = reflectionCoefficient(1,:);
        % %reflectionCoefficient = diag(reflectionCoefficient) * Jones;
        % %reflectionLoss = prod(abs(reflectionCoefficient).^2); 
        % reflectionLoss = -20*log10(reflectionLoss);
        E_inc1 = [0,1,0];
        E_inc2 = [1,0,0];
        for i = 1:orderReflection
            %is_perpendicular = norm(cross(E_inc, inc_plane(i)));
            is_parallel1 = abs(dot(E_inc1, inc_plane(i,:)));
            is_parallel2 = abs(dot(E_inc2, inc_plane(i,:)));
            if is_parallel1 < 1e-10 || is_parallel2< 1e-10
                R = [R,reflectionCoefficient(1,i)];
            else
                R = [R,reflectionCoefficient(2,i)];
            end
        end
        reflectionCoefficient = prod(sqrt(abs(R).^2));
        reflectionLoss = -20*log10(reflectionCoefficient);
    elseif strcmp(polarization,'dual')
        %strcmp(polarization,'dual') || strcmp(polarization,'cross')
        %reflectionCoefficient = prod(sqrt(abs(reflectionCoefficient).^2),2);
        % R = diag(reflectionCoefficient) * Jones;
        % reflectionLoss_H = -20*log10(R(1,1));
        % reflectionLoss_V = -20*log10(R(2,2));
        % reflectionLoss = [reflectionLoss_H;reflectionLoss_V];
        R = prod(sqrt(abs(reflectionCoefficient).^2),2);
        reflectionLoss = [-20*log10(R(1));-20*log10(R(2))];
    end
    diffuseLoss = [];
    scatterPaths = [];
    end
else
    %rou = exp(-0.5*(4*pi*r*cosd(aor)/lab)^2)*besselj(0,8*(pi*r*cos(aor)/lab));% Beckmann-Kirchhoff theory in the reflection direction
    scatteringcoefficient = 1 - Rou;
    %reflectionCoefficient = rou.*reflectionCoefficient;
    reflectionCoefficient = prod(sqrt(abs(reflectionCoefficient).^2),2);
    if strcmp(polarization,'V-V') || strcmp(polarization,'H-H')
        R = diag(reflectionCoefficient)*Jones;
        reflectionLoss = sum(R(:));
        %reflectionCoefficient = reflectionCoefficient(1,:);
        %reflectionCoefficient = diag(reflectionCoefficient) * Jones;
        %reflectionLoss = prod(abs(reflectionCoefficient).^2); 
        reflectionLoss_dB = -20 * log10(reflectionLoss);
        %diffuseLoss = -20 * log10(reflectionLoss / (prod(1-Rou) * (1-rou));
        diffuseLoss = reflectionLoss / prod(Rou)^2 * (1 - prod(Rou)^2);
        diffuseLoss_dB = -20 * log10(diffuseLoss);
        %directions = demoDiffuseScattering(30);
        %rayPower = -10 * log10(( 10 .^ (-diffuseLoss / 10) / 30) .* ones(30,1));
    elseif strcmp(polarization,'dual')
        %strcmp(polarization,'dual') || strcmp(polarization,'cross')
        %reflectionCoefficient = prod(sqrt(abs(reflectionCoefficient).^2),2);
        R = diag(reflectionCoefficient)*Jones;
        reflectionLoss_H_dB = -20*log10(R(1,1));
        diffuseLoss_H = R(1,1) / prod(Rou)^2 * (1 - prod(Rou)^2);
        diffuseLoss_H_dB = -20 * log10(diffuseLoss_H);
        %diffuseLoss_H = -20 * log10(R(1,1) / rou * (1-rou));
        reflectionLoss_V_dB = -20*log10(R(2,2));
        diffuseLoss_V = R(2,2) / prod(Rou)^2 * (1 - prod(Rou)^2);
        diffuseLoss_V_dB = -20 * log10(diffuseLoss_V);
        %diffuseLoss_V = -20 * log10(R(2,2) / rou * (1-rou));
        reflectionLoss = [reflectionLoss_H_dB;reflectionLoss_V_dB];
        diffuseLoss = [diffuseLoss_H_dB;diffuseLoss_V_dB];
        %directions = demoDiffuseScattering(30);
        %rayPower = -10 * log10(( 10 .^ (-diffuseLoss / 10) / 30) .* ones(2,30));
    end
    scatterPaths = computeDiffuseScattering(diffuseLoss,normal,doa,dod, 1 - prod(Rou)^2, randi([10,30]));
    % lx=0.01;
    % ly=0.01;
    % L=0.03;%correlation length
    % A=lx*ly;%illuminated area
    % if orderReflection==1
    %     theta1=Incident_angle;
    % else
    %     theta1=Incident_angle(2);
    % end
    % theta3=0;
    % k=2*pi/(3e8/140e9);
    % rou_star=[];
    % %for theta2=0:90
    % tmp=0:1:90;
    % thetaidx=[tmp Incident_angle];
    % for i=1:92
    %     theta2=thetaidx(i);
    %     g=k^2*r^2*(cosd(theta1)^2+cosd(theta2)^2);
    %     vx=k*(sind(theta1)-sind(theta2))*cosd(theta3);
    %     vy=k*(-sind(theta2)*sin(theta3));
    %     rou0=sinc(vx*lx)*sinc(vy*ly);
    %     F=(1+cosd(theta1)*cosd(theta2)-sind(theta1)*sind(theta2)*cosd(theta3));
    %     s=0;
    %     for i=1:1000
    %         s=s+g^i/(factorial(i)*i)*exp(-(vx^2+vy^2)*L^2/(4*i));
    %     end
    %     rou_inf=exp(-g)*(rou0^2+pi*L^2*F^2/A*s);
    %     if strcmp(polarization,'V-V')
    %         reflectionCoefficient1 = reflectionCoefficient(1,:)./rou;
    %         r_star=prod(abs(reflectionCoefficient1).^2); 
    %         rou_star=[rou_star r_star*rou_inf];
    %         reflectionLoss = -20*log10(rou_star);
    %     elseif strcmp(polarization,'H-H')
    %         reflectionCoefficient1 = reflectionCoefficient(2,:)./rou;
    %         r_star=prod(abs(reflectionCoefficient1).^2); 
    %         rou_star=[rou_star r_star*rou_inf];
    %         reflectionLoss = -20*log10(rou_star);
    %     end
    % end 
    %idx=92;
    %reflectionLoss=rou_star(idx);
    %reflectionLoss = -10*log10(rou_star);
end
