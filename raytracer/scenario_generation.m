function s = scenario_generation(s,confg_dim,IndoorSwitch)
% This function imports the configuration dimension and generates scatterers in the floor map

% Inputs:
% confg_dim:if indoor, room dimension. if outdoor, reference dimension
% s: the floor map which has been converted to struct format
% IndoorSwitch: If it's indoor, the walls, the ceiling and the floor will
% be added.
% x_cen/y_cen/z_cen: the centroid of scatterers
% x_dep/y_dep/z_dep: the size of scatterers
% scatterer_num: scatterer number

% Outputs:
% s: the floor map

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
if IndoorSwitch==1 
    for i=1:length(s.amf.object)
        for j=1:length(s.amf.object{1,i}.mesh.vertices.vertex)
            tmp=s.amf.object{1,i}.mesh.vertices.vertex{1,j}.coordinates;
            tmp_t=s.amf.object{1,i}.metadata.Text;
            tmp_x=tmp.x.Text;
            tmp_y=tmp.y.Text;
            tmp_z=tmp.z.Text;
            % Assign the configured dimension to the initial indoor floor map
            if (contains(tmp_t,'Wall') || contains(tmp_t,'Floor') || contains(tmp_t,'Ceiling'))
                if ~(contains(tmp_t,'Window')&&contains(tmp_t,'Door'))
                    if (~strcmp(tmp_x,'0'))
                        tmp.x.Text=num2str(confg_dim(1)*1000);%unit conversion
                    end
                    if (~strcmp(tmp_y,'0'))
                        tmp.y.Text=num2str(confg_dim(2)*1000);
                    end
                    if (~strcmp(tmp_z,'0'))
                        tmp.z.Text=num2str(confg_dim(3)*1000);
                    end
                    s.amf.object{1,i}.mesh.volume.triangle{1,1}.v1.Text;
                end
            end
            s.amf.object{1,i}.mesh.vertices.vertex{1,j}.coordinates=tmp;
        end
    end
else
    % for outdoor scenarios, clear the initial indoor floor map and add the
    % ground
    s.amf.object={};
    s.amf.material={};
    s.amf.object{1,1}.Attributes.id='0';
    s.amf.object{1,1}.metadata.Attributes.type='name';
    s.amf.object{1,1}.metadata.Text='Floor';
    s.amf.object{1,1}.mesh.vertices.vertex{1,1}.coordinates.x.Text=num2str(confg_dim(1)*1000);
    s.amf.object{1,1}.mesh.vertices.vertex{1,1}.coordinates.y.Text=num2str(confg_dim(2)*1000);
    s.amf.object{1,1}.mesh.vertices.vertex{1,1}.coordinates.z.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,4}.coordinates.x.Text=num2str(confg_dim(1)*1000);
    s.amf.object{1,1}.mesh.vertices.vertex{1,4}.coordinates.y.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,4}.coordinates.z.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,3}.coordinates.x.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,3}.coordinates.y.Text=num2str(confg_dim(2)*1000);
    s.amf.object{1,1}.mesh.vertices.vertex{1,3}.coordinates.z.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,2}.coordinates.x.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,2}.coordinates.y.Text=num2str(0);
    s.amf.object{1,1}.mesh.vertices.vertex{1,2}.coordinates.z.Text=num2str(0);
    s.amf.object{1,1}.mesh.volume.Attributes.materialid='1';
    s.amf.object{1,1}.mesh.volume.triangle{1,1}.v1.Text='0';
    s.amf.object{1,1}.mesh.volume.triangle{1,1}.v2.Text='1';
    s.amf.object{1,1}.mesh.volume.triangle{1,1}.v3.Text='2';
    s.amf.object{1,1}.mesh.volume.triangle{1,2}.v1.Text='3';
    s.amf.object{1,1}.mesh.volume.triangle{1,2}.v2.Text='1';
    s.amf.object{1,1}.mesh.volume.triangle{1,2}.v3.Text='0';
    s.amf.material{1,1}.Attributes.id='1';
    s.amf.material{1,1}.color.r.Text='0.1';
    s.amf.material{1,1}.color.g.Text='0.1';
    s.amf.material{1,1}.color.b.Text='0.1';
    s.amf.material{1,1}.metadata.Text='Floor';
end
end