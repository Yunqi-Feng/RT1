% generate the UE grid points:

function generate_UEgrid(num_RX, is_plot)

global s

if ~exist('is_plot','var')
     is_plot = 0;
end

M  = num_RX; % the number of UEs

Mx = s.Topology.UE_gridNum(1); % number of grid points alog x
My = s.Topology.UE_gridNum(2); % number of grid points alog y

if s.Topology.confg_dim(1) > 0
    xWidth = (s.Topology.confg_dim(1)-2*s.Topology.UE_safty(1))/Mx;
    if xWidth<2*s.Topology.UE_safty(1)
        error('inappropriate UE grid configuration along x-axis')
    end
    UE_x = s.Topology.UE_safty(1) + xWidth/2 + (0:Mx-1)*xWidth;
else
    Mx = 1;
    xWidth = 0;
    UE_x   = 0;
end
if s.Topology.confg_dim(2) > 0
    yWidth = (s.Topology.confg_dim(2)-2*s.Topology.UE_safty(2))/My;
    if yWidth<2*s.Topology.UE_safty(2)
        error('inappropriate UE grid configuration along y-axis')
    end
    UE_y = s.Topology.UE_safty(2) + yWidth/2 + (0:My-1)*yWidth;
else
    My     = 1;
    yWidth = 0;
    UE_y   = 0;
end
% xWidth = (s.Topology.confg_dim(1)-2*s.Topology.UE_safty(1))/Mx;
% yWidth = (s.Topology.confg_dim(2)-2*s.Topology.UE_safty(2))/My;
% 
% if xWidth<2*s.Topology.UE_safty(1)|| yWidth<2*s.Topology.UE_safty(2)
%     error('inappropriate UE grid configuration')
% end
% UE_x = s.Topology.UE_safty(1) + xWidth/2 + (0:Mx-1)*xWidth;
% UE_y = s.Topology.UE_safty(2) + yWidth/2 + (0:My-1)*yWidth;

if strcmp(s.Topology.UE_arch,'grid') || strcmp(s.Topology.UE_arch,'random')    
    if Mx*My < M
        error('The number of of grid points is smaller than the number of UEs...')
    elseif strcmp(s.Topology.UE_arch,'random') && (Mx*My == M)
        error('For random grid, you need more grid points  than the number of UEs...')
    elseif strcmp(s.Topology.UE_arch,'grid') && (Mx*My ~= M)
        error('For the (static) grid, you need the a number grid points  = the number of UEs...!')
    end
  
    newAllPoints = [];
    for ii=1:length(UE_x)
        for jj=1:length(UE_y)
            newAllPoints = [newAllPoints ; UE_x(ii) UE_y(jj)];
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
        s.Topology.UEgrid = newAllPoints;
        clear indices randomIndices points newAllPoints
    elseif length(newAllPoints)==M
        finalPoints = newAllPoints;
        clear newAllPoints
    else
        error('Number of grids is lower than that of APs');
    end
    rng('shuffle');
    
    if is_plot
        figure(650); movegui('northeast');
        subplot(2,1,1)
        if s.Topology.confg_dim(1)*s.Topology.confg_dim(2) > 0
            plot(s.Topology.UEgrid(:,1),s.Topology.UEgrid(:,2),'gO','MarkerSize',10,'MarkerFaceColor','g','DisplayName','UE grid');
            xlabel('x (m)'); ylabel('y (m)');  grid on; hold on; legend show
        elseif s.Topology.confg_dim(1) > 0
            plot(s.Topology.UEgrid(:,1),s.Topology.UE_height*ones(Mx*My,1),'gO','MarkerSize',10,'MarkerFaceColor','g','DisplayName','UE grid');
            xlabel('x (m)'); ylabel('z (m)');  grid on; hold on; legend show
        elseif s.Topology.confg_dim(2) > 0
            plot(s.Topology.UEgrid(:,2),s.Topology.UE_height*ones(Mx*My,1),'gO','MarkerSize',10,'MarkerFaceColor','g','DisplayName','UE grid');
            xlabel('y (m)'); ylabel('z (m)');  grid on; hold on; legend show
        end
        %
        subplot(2,1,2)
        plot3(s.Topology.UEgrid(:,1),s.Topology.UEgrid(:,2),s.Topology.UE_height*ones(Mx*My,1),'gO','MarkerSize',10,'MarkerFaceColor','g','DisplayName','UE grid');
        xlabel('x (m)'); ylabel('y (m)');  zlabel('z (m)'); grid on; hold on; legend show
    end
else
    error('UE topology not reconized...')
end

s.Topology.UE_x = finalPoints(:,1);
s.Topology.UE_y = finalPoints(:,2);
s.Topology.UE_z = s.Topology.UE_height * ones(M,1);

clear xWidth yWidth UE_x UE_y finalPoints;

