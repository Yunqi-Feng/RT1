function [x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,scatterer_num] = scatterer_config(confg_dim,scatterer_num,...
    scatterer_dim,scatterer_dis,NodeLocTx,NodeLocRx,IndoorSwitch, ...
    scatterer_confg,scatterer_height_max)
% This function generate scatterer properties based on number,
% distribution, dimension and location.

% Inputs:
% confg_dim:room dimension
% scatterer_num: scatterer dimension
% scatterer_dis: the statistical scatterer location distribution
% either (uniform/hyperbolic/Gaussian/Rayleigh)+(disk/sphere), or clusters
% NodeLocTx: a certain AP location
% NodeLocRx: a certain UE location
% IndoorSwitch: indoor/outdoor

% Outputs:
% x_cen/y_cen/z_cen: the centroid of scatterers
% x_dep/y_dep/z_dep: the size of scatterers

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% December 2023
x1=NodeLocTx(1);
y1=NodeLocTx(2);
z1=NodeLocTx(3);
x2=NodeLocRx(1);
y2=NodeLocRx(2);
z2=NodeLocRx(3);
x3=max(x1,x2);
y3=max(y1,y2);
z3=max(z1,z2);
x4=min(x1,x2);
y4=min(y1,y2);
z4=min(z1,z2);
d=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
x_cen=[];
y_cen=[];
z_cen=[];
switch(scatterer_dim)
    case 0
        x_trsh=min([0.5,0.5*confg_dim(1)]);
        y_trsh=min([0.5,0.5*confg_dim(2)]);
        z_trsh=min([0.25,0.7*confg_dim(3)]);
        %z_trsh=0.7*confg_dim(3);
    case 1
        x_trsh=min([2,0.5*confg_dim(1)]);
        y_trsh=min([2,0.5*confg_dim(2)]);
        z_trsh=min([1,0.7*confg_dim(3)]);
        %z_trsh=0.7*confg_dim(3);
    case 2
        x_trsh=20;
        y_trsh=20;
        %z_trsh=min([2.5,0.7*confg_dim(3)]);
        z_trsh=50;
        % x_trsh=min([15,0.5*confg_dim(1)]);
        % y_trsh=min([15,0.5*confg_dim(2)]);
        % z_trsh=20;
    case 3
        x_trsh=50;
        y_trsh=50;
        z_trsh=50;
    otherwise
        error("invalid scatterer dimension")
end
for i=1:scatterer_num
    if scatterer_dim~=3
        x_dep(i)=x_trsh*(0.5+0.5*rand);
        y_dep(i)=y_trsh*(0.5+0.5*rand);
        z_dep(i)=5+45*rand;
    else
        x_dep(i)=x_trsh*rand;
        y_dep(i)=y_trsh*rand;
        z_dep(i)=5+45*rand;
    end
end
int_safety=[max(x_dep) max(y_dep) max(z_dep)];
switch(scatterer_dim)
    case 0
        confg_safety=[0.5 0.5 0.25];
    case 1
        confg_safety=[1 1 0.5];
    case 2
        confg_safety=[5 5 0];
    case 3
        confg_safety=[5 5 0];
    otherwise
        error('invalid configuration')
end
%z_trsh1=adjust_trsh(scatterer_height_max,z_trsh,confg_safety,int_safety);
% if z_trsh1~=z_trsh
%     for i=1:scatterer_num
%         if scatterer_dim~=3
%             z_dep(i)=5+45*rand;
%         else
%             z_dep(i)=5+45*rand;
%         end
%     end
% end
% initialization
if scatterer_num>=1
    if isempty(scatterer_confg)
        for i=1:scatterer_num
            if scatterer_dis(1)==0%disk (ring/circle)
                dis=scatterer_dis(2);
                in=scatterer_dis(3);%if in=0, it's a circle; if in>0, it's a ring.
                ou=scatterer_dis(4);
                mean=scatterer_dis(5);
                std=scatterer_dis(6);
                [x_cen(i) y_cen(i) z_cen(i)]=disk(dis,in,ou,d,mean,std,NodeLocTx,NodeLocRx,scatterer_height_max);
                [x_cen(i),y_cen(i),z_cen(i)]=adjust_coordinates(x_cen(i),y_cen(i),z_cen(i),x_trsh,y_trsh,z_trsh,confg_dim,IndoorSwitch,scatterer_height_max);
            elseif scatterer_dis(1)==1%sphere
                dis=scatterer_dis(2);
                in=scatterer_dis(3);
                ou=scatterer_dis(4);
                mean=scatterer_dis(5);
                std=scatterer_dis(6);
                [x_cen(i),y_cen(i),z_cen(i)]=mysphere(dis,in,ou,d,mean,std,NodeLocTx,NodeLocRx,scatterer_height_max);
                [x_cen(i),y_cen(i),z_cen(i)]=adjust_coordinates(x_cen(i),y_cen(i),z_cen(i),x_dep,y_dep,z_dep,confg_dim,IndoorSwitch,scatterer_height_max);
            end
        end
        if scatterer_dis(1)==2%clusters
            cluster_num=scatterer_dis(2);%number of clusters
            scatterer_max=scatterer_dis(3);%maximum number of scatterers per cluster
            cluster_rad_max=scatterer_dis(4);
            clu_dis=scatterer_dis(5);
            sca_dis=scatterer_dis(6);
            [x_cen,y_cen,z_cen]=mycluster(cluster_num,scatterer_max,cluster_rad_max,clu_dis,sca_dis,scatterer_num,confg_dim,scatterer_confg,scatterer_height_max);
            for i=1:scatterer_num
                [x_cen(i),y_cen(i),z_cen(i)]=adjust_coordinates(x_cen(i),y_cen(i),z_cen(i),x_trsh,y_trsh,z_trsh,confg_dim,IndoorSwitch,scatterer_height_max);
            end     
        end
        if scatterer_dis(1)==0
            [x_cen,y_cen,z_cen]=generate_scatterergrid_2d(confg_dim,scatterer_num,scatterer_dim,IndoorSwitch,x_dep,y_dep,z_dep,x_cen,y_cen,z_cen,NodeLocTx,NodeLocRx,scatterer_height_max,int_safety,confg_safety);
        elseif scatterer_dis(1)==1 || scatterer_dis(1)==2
            [x_cen,y_cen,z_cen]=generate_scatterergrid_3d(confg_dim,scatterer_num,scatterer_dim,IndoorSwitch,x_dep,y_dep,z_dep,x_cen,y_cen,z_cen,scatterer_height_max,int_safety,confg_safety);
        elseif scatterer_dis(1)==3
            [x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,scatterer_num]=return_scatterergrid_3d(confg_dim,scatterer_num,scatterer_dim,x_dep,y_dep,z_dep,scatterer_height_max,IndoorSwitch);
            x_dep=x_dep';
            y_dep=y_dep';
            z_dep=z_dep';
        end
    elseif ~isempty(scatterer_confg)
        x_cen=scatterer_confg(1,:);
        y_cen=scatterer_confg(2,:);
        z_cen=scatterer_confg(3,:);
        x_dep=scatterer_confg(4,:);
        y_dep=scatterer_confg(5,:);
        z_dep=scatterer_confg(6,:);
    end
else
    x_cen=[];
    y_cen=[];
    z_cen=[];
    x_dep=[];
    y_dep=[];
    z_dep=[];
end
end