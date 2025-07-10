function [x_cen,y_cen,z_cen] = adjust_coordinates(x_cen,y_cen,z_cen,x_trsh, ...
    y_trsh,z_trsh,confg_dim,IndoorSwitch,scatterer_height_max)
% This function adjusts the the locations of scatterers to be practical
% Input:
% x_cen/y_cen/z_cen: the initial centroid of scatterers
% x_dep/y_dep/z_dep: the initial size of scatterers
% x_trsh/y_trsh/z_trsh: thresholds for the sizes
% confg_dim: the configured dimension
% IndoorSwitch: indoor/outdoor

%Output:
% x_cen/y_cen/z_cen: the modified centroid of scatterers
% x_dep/y_dep/z_dep: the modified size of scatterers

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
height_max=min([scatterer_height_max,confg_dim(3)-0.5*z_trsh]);
if (x_cen<0.5*x_trsh && IndoorSwitch==1)% indoor bottom
    x_cen=0.5*x_trsh;
elseif (x_cen>confg_dim(1)-0.5*x_trsh && IndoorSwitch==1)% indoor top
    x_cen=confg_dim(1)-0.5*x_trsh;
end
if (y_cen<0.5*y_trsh && IndoorSwitch==1)% indoor bottom
    y_cen=0.5*y_trsh;
elseif (y_cen>confg_dim(2)-0.5*y_trsh && IndoorSwitch==1)% indoor top
    y_cen=confg_dim(2)-0.5*y_trsh;
end
if z_cen<0.5*z_trsh% indoor/outdoor bottom
    z_cen=0.5*z_trsh;
elseif (z_cen>height_max)% indoor top
    z_cen=height_max;
end
end