function [x_cen,y_cen,z_cen] = mycluster(cluster_num,scatterer_max,cluster_rad_max,clu_dis, ...
    sca_dis,scatterer_num,confg_dim,scatterer_confg,scatterer_height_max)
if isempty(scatterer_confg)
    if (cluster_num>scatterer_num)
        error('inapproipriate configuration, cluster_num should be no less than scatterer_num')
    end
    if (scatterer_max<ceil(scatterer_num/cluster_num))
        error('inapproipriate configuration, please increase scatterer_max so that scatterer_max larger than scatterer_num/cluster_num')
    end

    if scatterer_num<3*cluster_num%At least 3 scatterers in a cluster
        error('scatterers are too few, at least 3 scatterers in a cluster')
    end
    if cluster_num==1
        scatterer_cluster=scatterer_num;
        cluster_rad_max=min([0.5*(confg_dim(3)) cluster_rad_max]);
        X=[];
        Y=[];
        Z=[];
        for i=1:length(scatterer_cluster)
            if clu_dis==0
                x=confg_dim(1)*rand;
                y=confg_dim(2)*rand;
                z=min(confg_dim(3)*rand,scatterer_height_max);
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==1
                std=1;
                bet=(1.5323/std^2)^(1/3);
                myrand=bet/2.367*sech((bet)^2*rand^2);
                x=confg_dim(1)*myrand;
                y=confg_dim(2)*myrand;
                z=min(confg_dim(3)*myrand,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=confg_dim(1)*myrand;
                end
                while y<0 || y>confg_dim(2)
                    y=confg_dim(2)*myrand;
                end
                while z<0 || z>confg_dim(3)
                    z=min(confg_dim(3)*myrand,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==2
                std=5;
                x=0.5*confg_dim(1)+std*randn;
                y=0.5*confg_dim(2)+std*randn;
                z=min(0.5*confg_dim(3)+std*randn,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=0.5*confg_dim(1)+std*randn;
                end
                while y<0 || y>confg_dim(2)
                    y=0.5*confg_dim(2)+std*randn;
                end
                while z<0 || z>confg_dim(3)
                    z=min(0.5*confg_dim(3)+std*randn,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==3
                tmp=rand;
                std=1;
                myrand=raylpdf(tmp,std);
                x=confg_dim(1)*myrand;
                y=confg_dim(2)*myrand;
                z=min(confg_dim(3)*myrand,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=confg_dim(1)*myrand;
                end
                while y<0 || y>confg_dim(2)
                    y=confg_dim(2)*myrand;
                end
                while z<0 || z>confg_dim(3)
                    z=min(confg_dim(3)*myrand,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            end
        end
        [X,Y,Z]=generate_clustergrid_3d(confg_dim,cluster_num,cluster_rad_max,X,Y,Z);
        Cluster_cor=[X;Y;Z]';
        x_cen=[];
        y_cen=[];
        z_cen=[];
        for i=1:length(scatterer_cluster)
            for j=1:scatterer_cluster(i)
                if sca_dis==0
                    r=cluster_rad_max*rand;
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=min(Cluster_cor(i,3)+r*cos(theta),scatterer_height_max);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==1
                    std=1;
                    bet=(1.5323/std^2)^(1/3);
                    myrand=bet/2.367*sech((bet)^2*rand^2);
                     r=cluster_rad_max*myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=bet/2.367*sech((bet)^2*(rand)^2);
                        r=cluster_rad_max*myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=min(Cluster_cor(i,3)+r*cos(theta),scatterer_height_max);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==2
                    std=5;
                    myrand=0.5*cluster_rad_max+std*randn;
                    r=myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=0.5*cluster_rad_max+std*randn;
                        r=myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=min(Cluster_cor(i,3)+r*cos(theta),scatterer_height_max);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==3
                    tmp=rand;
                    std=1;
                    myrand=raylpdf(tmp,std);
                    r=cluster_rad_max*myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=raylpdf(tmp,std);
                        r=cluster_rad_max*myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=min(Cluster_cor(i,3)+r*cos(theta),scatterer_height_max);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                end
            end
        end
    end
    if cluster_num>1
        for i=1:cluster_num-1
            tmp=randi(scatterer_max); 
            if tmp<3
                tmp=randi(scatterer_max);
            end
            scatterer_cluster{i}=tmp; 
        end
        tmp=sum(cell2mat(scatterer_cluster));
        scatterer_cluster{cluster_num}=scatterer_num-tmp;
        while scatterer_cluster{cluster_num}<3
            scatterer_cluster{cluster_num}=0;
            for i=1:cluster_num-1
                scatterer_cluster{i}=randi(scatterer_max);  
            end
            tmp=sum(cell2mat(scatterer_cluster));
            scatterer_cluster{cluster_num}=scatterer_num-tmp;
        end
        scatterer_cluster=cell2mat(scatterer_cluster);
        cluster_rad_max=min([0.5*(confg_dim(3)) cluster_rad_max]);
        X=[];
        Y=[];
        Z=[];
        for i=1:length(scatterer_cluster)
            if clu_dis==0
                x=confg_dim(1)*rand;
                y=confg_dim(2)*rand;
                z=min(confg_dim(3)*rand,scatterer_height_max);
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==1
                std=1;
                bet=(1.5323/std^2)^(1/3);
                myrand=bet/2.367*sech((bet)^2*rand^2);
                x=confg_dim(1)*myrand;
                y=confg_dim(2)*myrand;
                z=min(confg_dim(3)*myrand,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=confg_dim(1)*myrand;
                end
                while y<0 || y>confg_dim(2)
                    y=confg_dim(2)*myrand;
                end
                while z<0 || z>confg_dim(3)
                    z=min(confg_dim(3)*myrand,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==2
                std=5;
                x=0.5*confg_dim(1)+std*randn;
                y=0.5*confg_dim(2)+std*randn;
                z=min(0.5*confg_dim(3)+std*randn,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=0.5*confg_dim(1)+std*randn;
                end
                while y<0 || y>confg_dim(2)
                    y=0.5*confg_dim(2)+std*randn;
                end
                while z<0 || z>confg_dim(3)
                    z=min(0.5*confg_dim(1)+std*randn,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            elseif clu_dis==3
                tmp=rand;
                std=1;
                myrand=raylpdf(tmp,std);
                x=confg_dim(1)*myrand;
                y=confg_dim(2)*myrand;
                z=min(confg_dim(3)*myrand,scatterer_height_max);
                while x<0 || x>confg_dim(1)
                    x=confg_dim(1)*myrand;
                end
                while y<0 || y>confg_dim(2)
                    y=confg_dim(2)*myrand;
                end
                while z<0 || z>confg_dim(3)
                    z=min(confg_dim(3)*myrand,scatterer_height_max);
                end
                X=[X x];
                Y=[Y y];
                Z=[Z z];
            end
        end
        [X,Y,Z]=generate_clustergrid_3d(confg_dim,cluster_num,cluster_rad_max,X,Y,Z);
        Cluster_cor=[X;Y;Z]';
        x_cen=[];
        y_cen=[];
        z_cen=[];
        for i=1:length(scatterer_cluster)
            for j=1:scatterer_cluster(i)
                if sca_dis==0
                    r=cluster_rad_max*rand;
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=Cluster_cor(i,3)+r*cos(theta);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==1
                    std=1;
                    bet=(1.5323/std^2)^(1/3);
                    myrand=bet/2.367*sech((bet)^2*rand^2);
                     r=cluster_rad_max*myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=bet/2.367*sech((bet)^2*(rand)^2);
                        r=cluster_rad_max*myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=Cluster_cor(i,3)+r*cos(theta);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==2
                    std=5;
                    myrand=0.5*cluster_rad_max+std*randn;
                    r=myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=0.5*cluster_rad_max+std*randn;
                        r=myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=Cluster_cor(i,3)+r*cos(theta);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                elseif sca_dis==3
                    tmp=rand;
                    std=1;
                    myrand=raylpdf(tmp,std);
                    r=cluster_rad_max*myrand;
                    while(r>cluster_rad_max)|| (r<0)
                        myrand=raylpdf(tmp,std);
                        r=cluster_rad_max*myrand;
                    end
                    phi=2*pi*rand;
                    theta=-pi+2*pi*rand;
                    x_cen_tmp=Cluster_cor(i,1)+r*sin(theta)*cos(phi);
                    y_cen_tmp=Cluster_cor(i,2)+r*sin(theta)*sin(phi);
                    z_cen_tmp=Cluster_cor(i,3)+r*cos(theta);
                    x_cen=[x_cen x_cen_tmp];
                    y_cen=[y_cen y_cen_tmp];
                    z_cen=[z_cen z_cen_tmp];
                end
            end
        end
    end
else
    x_cen=scatterer_confg(1,1:scatterer_num);
    y_cen=scatterer_confg(2,1:scatterer_num);
    z_cen=scatterer_confg(3,1:scatterer_num);
end
end