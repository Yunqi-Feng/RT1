% generate the AP grid points:

function  generate_APgrid(num_TX, is_plot)

global s

if ~exist('is_plot','var')
     is_plot = 0;
end

M = num_TX;% the number of APs

Mx = s.Topology.AP_gridNum(1); % number of grid points alog x
My = s.Topology.AP_gridNum(2); % number of grid points alog y

if ~strcmp(s.Topology.AP_arch,'ceil-collocated')
    if s.Topology.confg_dim(1) > 0
        xWidth = (s.Topology.confg_dim(1)-2*s.Topology.AP_safty(1))/Mx;
        if xWidth<2*s.Topology.AP_safty(1)
            error('inappropriate AP grid configuration along x-axis')
        end
        AP_x = s.Topology.AP_safty(1) + xWidth/2 + (0:Mx-1)*xWidth;
    else
        Mx     = 1;
        xWidth = 0;
        AP_x   = 0;
    end
    if s.Topology.confg_dim(2) > 0
        yWidth = (s.Topology.confg_dim(2)-2*s.Topology.AP_safty(2))/My;
        if yWidth<2*2*s.Topology.AP_safty(2)
            error('inappropriate AP grid configuration along y-axis')
        end
        AP_y = s.Topology.AP_safty(2) + yWidth/2 + (0:My-1)*yWidth;
    else
        My     = 1;
        yWidth = 0;
        AP_y   = 0;
    end
    %     if xWidth<2*s.Topology.AP_safty(1)|| yWidth<2*2*s.Topology.AP_safty(2)
    %         error('inappropriate AP grid configuration')
    %     end
    %     AP_x = s.Topology.AP_safty(1) + xWidth/2 + (0:Mx-1)*xWidth;
    %     AP_y = s.Topology.AP_safty(2) + yWidth/2 + (0:My-1)*yWidth;
end

if strcmp(s.Topology.AP_arch,'ceil-piazza')
    if 2*(Mx+My) ~= M
        error('The number of points in piazza is not equal to the number of APs!')
    end
    
    finalPoints = [ [AP_x'                  AP_y(1)*ones(Mx,1)]   ; ...
                     [AP_x(end)*ones(My,1)  AP_y']                ; ...
                     [AP_x(end:-1:1)'       AP_y(end)*ones(Mx,1)] ; ...
                     [AP_x(1)*ones(My,1)    AP_y(end:-1:1)']     ];
elseif strcmp(s.Topology.AP_arch,'ceil-grid') || strcmp(s.Topology.AP_arch,'ceil-random')    
    if Mx*My < M
        error('The number of of grid points is smaller than the number of APs...')
    elseif strcmp(s.Topology.AP_arch,'ceil-random') && (Mx*My == M)
        error('For random grid, you need more grid points  than the number of APs...')
    elseif strcmp(s.Topology.AP_arch,'ceil-grid') && (Mx*My ~= M)
        error('For the (static) grid, you need the a number grid points  must be = the number of APs...')
    end
  
    newAllPoints = [];
    for ii=1:length(AP_x)
        for jj=1:length(AP_y)
            newAllPoints = [newAllPoints ; AP_x(ii) AP_y(jj)];
        end
    end
    
     %s_rng = rng;
     rng('shuffle');
    if length(newAllPoints) > M %Random strategy
        [~,indices]       = sort(rand(size(newAllPoints,1),1));
        randomIndices=randperm(length(indices),M);
        %randomIndices = indices(1:M);
        points        = newAllPoints;
        finalPoints = points(randomIndices,:);
        s.Topology.APgrid = newAllPoints;
        clear indices randomIndices points newAllPoints
    elseif length(newAllPoints)==M
        finalPoints = newAllPoints;
        clear newAllPoints
    else
        error('Number of grid points is lower than the number of APs...');
    end
    rng('shuffle');
    
    if is_plot
        figure(650); clf; movegui('northeast');
        subplot(2,1,1)
        if s.Topology.confg_dim(1)*s.Topology.confg_dim(2) > 0
            plot(s.Topology.APgrid(:,1),s.Topology.APgrid(:,2),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('x (m)'); ylabel('y (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(1))]);
            ylim([0,max(1,s.Topology.confg_dim(2))]);
        elseif s.Topology.confg_dim(1) > 0
            plot(s.Topology.APgrid(:,1),s.Topology.AP_height*ones(Mx*My,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('x (m)'); ylabel('z (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(1))]);
            ylim([0,max(1,s.Topology.confg_dim(3))]);
        elseif s.Topology.confg_dim(2) > 0
            plot(s.Topology.APgrid(:,2),s.Topology.AP_height*ones(Mx*My,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('y (m)'); ylabel('z (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(2))]);
            ylim([0,max(1,s.Topology.confg_dim(3))]); 
        end
        subplot(2,1,2)
        plot3(s.Topology.APgrid(:,1),s.Topology.APgrid(:,2),s.Topology.AP_height*ones(Mx*My,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
        xlabel('x (m)'); ylabel('y (m)');  zlabel('z (m)'); grid on; hold on; 
        xlim([0,max(1,s.Topology.confg_dim(1))]); 
        ylim([0,max(1,s.Topology.confg_dim(2))]);
    end  
elseif strcmp(s.Topology.AP_arch,'ceil-collocated')
    s.Topology.AP_x = ones(1,M)  * s.Topology.confg_dim(1) / 2;
    s.Topology.AP_y = ones(1,M)  * s.Topology.confg_dim(2) / 2;
    
     if is_plot
        figure(650); clf; movegui('northeast'); 
        subplot(2,1,1)
        if s.Topology.confg_dim(1)*s.Topology.confg_dim(2) > 0
            plot(s.Topology.AP_x,s.Topology.AP_y,'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('x (m)'); ylabel('y (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(1))]);
            ylim([0,max(1,s.Topology.confg_dim(2))]);
        elseif s.Topology.confg_dim(1) > 0
            plot(s.Topology.AP_x,s.Topology.AP_height*ones(1*1,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('x (m)'); ylabel('z (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(1))]);
            ylim([0,max(1,s.Topology.confg_dim(3))])
        elseif s.Topology.confg_dim(2) > 0
            plot(s.Topology.AP_y,s.Topology.AP_height*ones(1*1,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
            xlabel('y (m)'); ylabel('z (m)');  grid on; hold on;
            xlim([0,max(1,s.Topology.confg_dim(2))]);
            ylim([0,max(1,s.Topology.confg_dim(3))])
        end
        subplot(2,1,2)
        plot3(s.Topology.AP_x,s.Topology.AP_y,s.Topology.AP_height*ones(1*1,1),'rO','MarkerSize',20,'MarkerFaceColor','r','DisplayName','AP grid');
        xlabel('x (m)'); ylabel('y (m)');  zlabel('z (m)'); grid on; hold on; 
        xlim([0,max(1,s.Topology.confg_dim(1))]); 
        ylim([0,max(1,s.Topology.confg_dim(2))]);
    end  
else
    error('AP topology not recofnized...')
end
if ~strcmp(s.Topology.AP_arch,'ceil-collocated')
    s.Topology.AP_x = finalPoints(:,1);
    s.Topology.AP_y = finalPoints(:,2);
end
clear xWidth yWidth AP_x AP_y finalPoints;

