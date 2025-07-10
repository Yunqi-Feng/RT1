function  [x_cen,y_cen,z_cen,x_dep,y_dep,z_dep,scatterer_num]=return_scatterergrid_3d(confg_dim,scatterer_num,scatterer_dim, ...
    x_dep,y_dep,z_dep,scatterer_height_max,IndoorSwitch)
% This function sets scatterer grid over the configured space. Safety
% distances are set for different scatterer dimensions. The Eventual grid
% points (the center of the scatterer) are selected from the points that
% are closest to initial scatterers.

% Input:
% confg_dim: the indoor/outdoor dimension we want to configure
% scatterer_num: scatterer number
% scatterer_dim: scatterer dimension, 0=small, 1=medium, 2=big
% IndoorSwitch: 0=outdoor, 1=outdoor
% x_cen, y_cen, z_cen: center of scatterers
% x_dep, y_dep, z_dep: dimension of scatterers

% Output:
% x_cen, y_cen, z_cen: modified centers of scatterers

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
for iii=1:1
    int_safety=[max(x_dep) max(y_dep) max(z_dep)];
    int_safety=[sqrt(int_safety(1)^2+int_safety(2)^2),sqrt(int_safety(1)^2+int_safety(2)^2),int_safety(3)];
    switch(scatterer_dim)
        case 0
            confg_safety=[0.5 0.5 0.25];
        case 1
            confg_safety=[1 1 0.5];
        case 2
            confg_safety=[2 2 0];
        case 3
            confg_safety=[2 2 0];
        otherwise
            error('invalid configuration')
    end
    %scatterer_safety=confg_safety+int_safety;
    scatterer_safety = confg_safety + int_safety;
    if IndoorSwitch==1
        if confg_dim(3)-scatterer_safety(3)<0 || confg_dim(2)-scatterer_safety(2)<0 || confg_dim(1)-scatterer_safety(1)<0 || confg_dim(3)-0.5*scatterer_safety(3)<0
            error('Please enalrge the room dimension.')
        end
        height_max=min([scatterer_height_max,confg_dim(3)-0.5*scatterer_safety(3)]);
    else
        if confg_dim(2)-scatterer_safety(2)<0 || confg_dim(1)-scatterer_safety(1)<0
            error('Please enalrge the configuration dimension.')
        end
        height_max=scatterer_height_max;
    end
    scatterer_gridNum=[floor(confg_dim(1)/(scatterer_safety(1))) floor(confg_dim(2)/(scatterer_safety(2))) floor((height_max+0.5*scatterer_safety(3))/(scatterer_safety(3)))];
    %scatterer_gridNum=[floor(confg_dim(1)/(1.25*max(x_dep))) floor(confg_dim(2)/(1.25*max(y_dep))) floor((height_max+0.5*scatterer_safety(3))/(scatterer_safety(3)))];
    Mx=scatterer_gridNum(1);
    My=scatterer_gridNum(2);
    %Mz=1;
    Mz=scatterer_gridNum(3);
    xWidth=scatterer_safety(1);
    yWidth=scatterer_safety(2);
    zWidth=scatterer_safety(3);  
    %xWidth=max(x_dep);
    %yWidth=max(y_dep);
    %zWidth=max(z_dep);
    scatterer_x=scatterer_safety(1)/2+(0:Mx-1)*xWidth;
    scatterer_y=scatterer_safety(2)/2+(0:My-1)*yWidth;
    scatterer_z=scatterer_safety(3)/2+(0:Mz-1)*zWidth; 
    newAllPoints = [];
    for ii=1:length(scatterer_x)
        for jj=1:length(scatterer_y)
            for kk=1:length(scatterer_z)
                newAllPoints = [newAllPoints;scatterer_x(ii) scatterer_y(jj) scatterer_z(kk)];
            end
        end
    end
    x_cen=newAllPoints(:,1);
    y_cen=newAllPoints(:,2);
    z_cen=newAllPoints(:,3);
    if scatterer_num>length(x_cen)
        error('The current maxixmum posible number of scatterers is %d, please reduce the configured scatterer number',length(x_cen))
    end
    [~,indices]   = sort(rand(size(newAllPoints,1),1));
    randomIndices = indices(1:scatterer_num);
    points=newAllPoints;
    finalPoints=points(randomIndices,:);
    %scatterer_num=length(x_cen);
    x_cen=finalPoints(:,1);
    y_cen=finalPoints(:,2);
    z_cen=finalPoints(:,3);
    if scatterer_dim == 2
        for i=1:scatterer_num
            z_cen(i)=0.5*z_dep(i);
        end
    end

    % for i=1:scatterer_num
    %     if scatterer_dim~=3
    %         x_dep(i)=x_trsh*(0.5+0.5*rand);
    %         y_dep(i)=y_trsh*(0.5+0.5*rand);
    %         z_dep(i)=z_trsh*(0.5+0.5*rand);
    %     else
    %         x_dep(i)=x_trsh*rand;
    %         y_dep(i)=y_trsh*rand;
    %         z_dep(i)=z_trsh*rand;
    %     end
    %end
end
end