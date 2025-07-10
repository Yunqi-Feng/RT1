function  [x_cen,y_cen,z_cen]=generate_scatterergrid(confg_dim,scatterer_num,scatterer_dim,scatterer_dis,IndoorSwitch,x_dep,y_dep,z_dep,x_cen,y_cen,z_cen,NodeLocRx)
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
if IndoorSwitch==1
    int_safety=[max(x_dep) max(y_dep) max(z_dep)];
    switch(scatterer_dim)
        case 0
            confg_safety=[0.5 0.5 0.25];
        case 1
            confg_safety=[1 1 0.5];
        case 2
            confg_safety=[2 2 1];
        case 3
            confg_safety=[2 2 1];
        otherwise
            error('invalid configuration')
    end
    scatterer_safety=confg_safety+int_safety;
    scatterer_gridNum=[floor((confg_dim(1)-scatterer_safety(1))/scatterer_safety(1))+1 floor((confg_dim(2)-scatterer_safety(2))/scatterer_safety(2))+1 floor((confg_dim(3)-scatterer_safety(3))/scatterer_safety(3))+1];
    Mx=scatterer_gridNum(1);
    My=scatterer_gridNum(2);
    Mz=scatterer_gridNum(3);
    if scatterer_num>Mx*My*Mz
        error("inappropriate configuration, too many scatterers.")
    end
    xWidth=(confg_dim(1)-scatterer_safety(1))/Mx;
    yWidth=(confg_dim(2)-scatterer_safety(2))/My;
    zWidth=(confg_dim(3)-scatterer_safety(3))/Mz;
    scatterer_x=scatterer_safety(1)/2+xWidth/2+(0:Mx-1)*xWidth;
    scatterer_y=scatterer_safety(2)/2+yWidth/2+(0:My-1)*yWidth;
    scatterer_z=scatterer_safety(3)/2+zWidth/2+(0:Mz-1)*zWidth;  
    newAllPoints = [];
    for ii=1:length(scatterer_x)
        for jj=1:length(scatterer_y)
            for kk=1:length(scatterer_z)
                newAllPoints = [newAllPoints;scatterer_x(ii) scatterer_y(jj) scatterer_z(kk)];
            end
        end
    end                                                                                                                                                         
    prePoints=[x_cen;y_cen;z_cen]';
    for i=1:scatterer_num
        dist=[];
        for j=1:size(newAllPoints,1)
            dist=[dist norm(prePoints(i,:)-newAllPoints(j,:))];
        end
        idx=find(dist==min(dist));
        if ~isempty(idx)
            idx=idx(1);
            x_cen(i)=newAllPoints(idx,1);
            y_cen(i)=newAllPoints(idx,2);
            z_cen(i)=newAllPoints(idx,3);
            newAllPoints(idx,:)=[];
        end    
    end
else
    int_safety=[max(x_dep) max(y_dep) max(z_dep)];
    switch(scatterer_dim)
        case 0
            confg_safety=[0.5 0.5 0.25];
        case 1
            confg_safety=[1 1 0.5];
        case 2
            confg_safety=[2 2 1];
        case 3
            confg_safety=[2 2 1];
        otherwise
            error("invalid configuration")
    end
    scatterer_safety=confg_safety+int_safety;
    dim1_max=max([max(x_cen)-min(x_cen) max(x_cen) confg_dim(1)]);
    dim2_max=max([max(y_cen)-min(y_cen) max(y_cen) confg_dim(2)]);
    dim3_max=max([max(z_cen)-min(z_cen) max(z_cen) confg_dim(3)]);
    dim1_min=min([min(x_cen) 0]);
    dim2_min=min([min(y_cen) 0]);
    dim3_min=min([min(z_cen) 0]);
    confg_dim=[dim1_max-dim1_min dim2_max-dim2_min dim3_max-dim3_min];
    scatterer_gridNum=[floor((confg_dim(1)-scatterer_safety(1))/scatterer_safety(1))+1 floor((confg_dim(2)-scatterer_safety(2))/scatterer_safety(2))+1 floor((confg_dim(3)-scatterer_safety(3))/scatterer_safety(3))+1];
    Mx=scatterer_gridNum(1);
    My=scatterer_gridNum(2);
    Mz=scatterer_gridNum(3);
    if scatterer_num>Mx*My*Mz
        error("inappropriate configuration, too many scatterers.")
    end
    xWidth=(confg_dim(1)-scatterer_safety(1))/Mx;
    yWidth=(confg_dim(2)-scatterer_safety(2))/My;
    zWidth=(confg_dim(3)-scatterer_safety(3))/Mz;
    scatterer_x=dim1_min+scatterer_safety(1)/2+xWidth/2+(0:Mx-1)*xWidth;
    scatterer_y=dim2_min+scatterer_safety(2)/2+yWidth/2+(0:My-1)*yWidth;
    scatterer_z=dim3_min+scatterer_safety(3)/2+zWidth/2+(0:Mz-1)*zWidth;  
    newAllPoints = [];
    for ii=1:length(scatterer_x)
        for jj=1:length(scatterer_y)
            for kk=1:length(scatterer_z)
                newAllPoints = [newAllPoints;scatterer_x(ii) scatterer_y(jj) scatterer_z(kk)];
            end
        end
    end
    %[~,indices]=sort(rand(size(newAllPoints,1),1));
    %randomIndices=indices(1:scatterer_num);
    %points=newAllPoints;
    %finalPoints=points(randomIndices,:);
    prePoints=[x_cen;y_cen;z_cen]';
    for i=1:scatterer_num
        dist=[];
        for j=1:size(newAllPoints,1)
            dist=[dist norm(prePoints(i,:)-newAllPoints(j,:))];
        end
        idx=find(dist==min(dist));
        if ~isempty(idx)
            x_cen(i)=newAllPoints(idx(1),1);
            y_cen(i)=newAllPoints(idx(1),2);
            z_cen(i)=newAllPoints(idx(1),3);
            newAllPoints(idx,:)=[];
        end
    end
end
end

