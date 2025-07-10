function  [x_cen,y_cen,z_cen]=generate_clustergrid_3d(confg_dim,cluster_num,cluster_rad_max,x_cen,y_cen,z_cen)
for ii=1:1
    cluster_safety=[cluster_rad_max cluster_rad_max cluster_rad_max];
    cluster_gridNum=[floor((confg_dim(1)-cluster_safety(1))/cluster_safety(1))+1 floor((confg_dim(2)-cluster_safety(2))/cluster_safety(2))+1 floor((confg_dim(3)-cluster_safety(3))/cluster_safety(3))+1];
    Mx=cluster_gridNum(1);
    My=cluster_gridNum(2);
    Mz=cluster_gridNum(3);
    if cluster_num>Mx*My*Mz
        error("inappropriate configuration, too many clusters.")
    end
    xWidth=(confg_dim(1)-cluster_safety(1))/Mx;
    yWidth=(confg_dim(2)-cluster_safety(2))/My;
    zWidth=(confg_dim(3)-cluster_safety(3))/Mz;
    cluster_x=cluster_safety(1)/2+xWidth/2+(0:Mx-1)*xWidth;
    cluster_y=cluster_safety(2)/2+yWidth/2+(0:My-1)*yWidth;
    cluster_z=cluster_safety(3)/2+zWidth/2+(0:Mz-1)*zWidth;  
    newAllPoints = [];
    for ii=1:length(cluster_x)
        for jj=1:length(cluster_y)
            for kk=1:length(cluster_z)
                newAllPoints = [newAllPoints;cluster_x(ii) cluster_y(jj) cluster_z(kk)];
            end
        end
    end                                                                                                                                                         
    prePoints=[x_cen;y_cen;z_cen]';
    for i=1:cluster_num
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
end
end

