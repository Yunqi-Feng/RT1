function [x_cen,y_cen,z_cen] = disk(dis,in,ou,d,mean,std,NodeLocTx, ...
    NodeLocRx,scatterer_height_max)
% This function generates scatterer coordinates for a statistical disk 
% distribution.

% Input: 
% dis:statistical distribution 
% in: inner radius
% ou: outer radius
% d: TX/RX didtance
% mean/std: for specific statistical distributions
% NodeLocTx: a certain AP location
% NodeLocRx: a certain UE location

% Outputs:
% x_cen/y_cen/z_cen: the centroid of scatterers
% x_dep/y_dep/z_dep: the size of scatterers

% Author: Yunqi Feng <Yunqi.Feng@UGent.be>
% September 2023
x1=NodeLocTx(1);
y1=NodeLocTx(2);
z1=NodeLocTx(3);
x2=NodeLocRx(1);
y2=NodeLocRx(2);
z2=NodeLocRx(3);
if mean<in || mean>ou
    error('inappropriate configuration, the mean of a disk should be within the inner radius and outer radius')
end
if in<0 || ou<0
    error('inappropriate configuration with negative radius')
%elseif d<ou && d>in
    %ou=d;
    %if mean>d
        %mean=0.5*(ou-in)+in;
    %end
%elseif in>d
    %in=d;
    %ou=2*d;
    %mean=0.5*(ou-in)+in;
end
if dis==0%uniform
    r=in+(ou-in)*rand;
    phi=2*pi*rand;
    x_cen=x2+r*cos(phi);
    y_cen=y2+r*sin(phi);
    z_cen=min(scatterer_height_max,z2);
elseif dis==1%hyperbolic
    bet=(1.5323/std^2)^(1/3);
    myrand=bet/2.367*sech((bet)^2*(rand)^2);
    r=in+(ou-in)*myrand;
    %r=myrand;
    while(r>ou)||(r<in)
        myrand=bet/2.367*sech((bet)^2*(rand)^2);
        r=in+(ou-in)*myrand;
    end
    phi=2*pi*rand;
    x_cen=x2+r*cos(phi);
    y_cen=y2+r*sin(phi);
    z_cen=min(scatterer_height_max,z2);
elseif dis==2%Gaussian
    myrand=mean+std*randn;
    r=myrand;
    %r=myrand;
    while(r>ou)||(r<in)
        myrand=mean+std*randn;
        r=myrand;
    end
    phi=2*pi*rand;
    x_cen=x2+r*cos(phi);
    y_cen=y2+r*sin(phi);
    z_cen=min(scatterer_height_max,z2);
elseif dis==3%Rayleigh
    tmp=rand;
    myrand=raylpdf(tmp,std);
    r=in+(ou-in)*myrand;
    while(r>ou)||(r<in)
        tmp=rand;
        myrand=raylpdf(tmp,std);
        r=in+(ou-in)*myrand;
    end
    phi=2*pi*rand;
    x_cen=x2+r*cos(phi);
    y_cen=y2+r*sin(phi);
    z_cen=min(scatterer_height_max,z2);
end