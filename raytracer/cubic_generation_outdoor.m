function [s] = cubic_generation_outdoor(s,idx,x_cen,x_dep,y_cen,y_dep,z_cen,z_dep)
% This function generates cubic scatterers according to the vertex
% coordinates for outdoor scenaros which have a ground.

% Input:
% s: the floor map which has been converted to struct format
% idx:loop index determined by numberof scatterers.
% x_cen, y_cen, z_cen: centers of scatterers
% x_dep, y_dep, z_dep: dimension of scatterers

% Output:
% s: the floor map

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
tmp=randi([0 1]);% control the rotation
if tmp==0
    %% scatterer vertex 0
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 1
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 2
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 3
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 4
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 5
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 6
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 7
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
elseif tmp==1
    ang=randi([0 90]);
    %% scatterer vertex 0
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2*cosd(ang)+y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2*cosd(ang)+x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,1}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 1
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2*cosd(ang)-y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2*cosd(ang)+x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,2}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 2
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2*cosd(ang)-y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2*cosd(ang)+x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,3}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 3
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.x.Text=num2str(1000*(x_cen(idx)+x_dep(idx)/2*cosd(ang)+y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2*cosd(ang)+x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,4}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 4
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2*cosd(ang)+y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2*cosd(ang)-x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,5}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 5
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2*cosd(ang)-y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2*cosd(ang)-x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,6}.coordinates.z.Text=num2str(1000*(z_cen(idx)-z_dep(idx)/2));
    %% scatterer vertex 6
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2*cosd(ang)-y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.y.Text=num2str(1000*(y_cen(idx)+y_dep(idx)/2*cosd(ang)-x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,7}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
    %% scatterer vertex 7
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.x.Text=num2str(1000*(x_cen(idx)-x_dep(idx)/2*cosd(ang)+y_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.y.Text=num2str(1000*(y_cen(idx)-y_dep(idx)/2*cosd(ang)-x_dep(idx)/2*sind(ang)));
    s.amf.object{1,1+idx}.mesh.vertices.vertex{1,8}.coordinates.z.Text=num2str(1000*(z_cen(idx)+z_dep(idx)/2));
end