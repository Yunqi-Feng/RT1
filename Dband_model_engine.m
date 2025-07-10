%% DBAND_MODEL_ENGINE generates channel impulse respone realizations based
% on an input configuration
%
%  Inputs:
%      config - configuration struct
%  Outputs:
%      output - cell array with CIR data including amplitude and delay
%      nodeLocTx - cell array with locations of transmit antennas
%      nodeLocRx - cell array with locations of receiving antennas
%
%  Author:  B. De Beelde <Brecht.DeBeelde@UGent.be>
%  Version: 1.2
%  September 2021; Last revision: June 17, 2022
%  Modified by Yunqi Feng <Yunqi.Feng@UGent.be> for a flexible AP/UE
%  configuration.
%  December 2023
function [output,scatterer_confg] = Dband_model_engine()
cd functions/Channel_model/DBARC/ % Q++: normally the addpath(genpath('functions')) should be enough ? 
%addpath(genpath('functions'))
global s
config = s.Dband_chConfig;
% Input System Parameters
if ~isfield(s,'paraCfg')
    s.paraCfg = getSystemParameters(config);     
end
paraCfg = s.paraCfg;
num=paraCfg.scatterer_num;
scatterer_confg=[];
% - Generate NEW AP location in case of random AP distribution
% AND Generate a randomly sorted UE grid points:
%   ==> For a single AP: we generate miltiple UE grid points that will be
%       used to pikup Kue points fulfilling the 'config.Topology.channel_req' 
if ~isfield(config.Topology,'nodeLocTx') 
    if strcmp(config.Topology.AP_arch,'ceil-random')
        newAllPoints  =  config.Topology.APgrid;
        [~,indices]   = sort(rand(size(newAllPoints,1),1));
        randomIndices = indices(1:config.num_TX);
        points        = newAllPoints;
        finalPoints   = points(randomIndices,:);
        s.chan_config.Topology.AP_x = finalPoints(:,1);
        s.chan_config.Topology.AP_y = finalPoints(:,2);
        clear indices randomIndices points newAllPoints finalPoints
    elseif strcmp(config.Topology.AP_arch,'ceil-collocated')
        s.chan_config.Topology.AP_x  = s.Topology.AP_x;
        s.chan_config.Topology.AP_y  = s.Topology.AP_y;
    end
    if strcmp(config.Topology.UE_arch,'random')
        newAllPoints =  config.Topology.UEgrid;
        [~,indices]  = sort(rand(size(newAllPoints,1),1));
        %s.chan_config.Topology.tmpUEgrid_idx=indices(1)
        s.chan_config.Topology.tmpUEgrid = newAllPoints(indices,:);
        clear newAllPoints indices
    else % grid case
        s.chan_config.Topology.tmpUEgrid = [s.Topology.UE_x s.Topology.UE_y];
    end
    s.chan_config.Topology.tmpUEgrid_idx = 1;
end
UE_points  = s.chan_config.Topology.tmpUEgrid;
%UEgrid_idx = s.chan_config.Topology.tmpUEgrid_idx;
Dz = config.Topology.confg_dim(3);

if config.is_plot 
    figure(651); movegui('northeast');  
end
UEgrid_idx         = s.chan_config.Topology.tmpUEgrid_idx;
Num_UE_grid_points = size(UE_points,1);
if config.num_RX > Num_UE_grid_points
    error('not enough UE grid points');
end
ref_nodeLocTx=[];
ref_nodeLocRx=[];
for ii=1:config.num_TX
    % Generate the AP locations: 
    if ~isfield(config.Topology,'nodeLocTx')
        x =  s.chan_config.Topology.AP_x(ii);
        y =  s.chan_config.Topology.AP_y(ii); 
        %x =  s.chan_config.Topology.AP_x(:,ii);
        %y =  s.chan_config.Topology.AP_y(:,ii); 
        if isfield(config.Topology,'AP_height')
            z = config.Topology.AP_height;
        else
            z = unifrnd(0.1,Dz);
        end
        %z=z*ones(config.Topology.timedivision,1);
        s.chan_config.Topology.nodeLocTx.x = x;
        s.chan_config.Topology.nodeLocTx.y = y;
        s.chan_config.Topology.nodeLocTx.z = z;
        s.Topology.AP_x=s.chan_config.Topology.AP_x;
        s.Topology.AP_y=s.chan_config.Topology.AP_y;
        % some intermediate plots:
        if config.is_plot
            figure(651);  
            subplot(2,1,1); plot(x,y,'O','MarkerSize',10,'MarkerFaceColor','r'); hold on;
            xlabel('x (m)'); ylabel('y (m)'); grid minor; xlim([0,s.Topology.confg_dim(1)]); ylim([0,s.Topology.confg_dim(2)]);
            subplot(2,1,2);  plot3(x,y,z,'O','MarkerSize',10,'MarkerFaceColor','r'); hold on;
            xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); grid minor; xlim([0,s.Topology.confg_dim(1)]); ylim([0,s.Topology.confg_dim(2)]);
        end
        clear x y z
    end
    nodeLocTx{ii}=s.chan_config.Topology.nodeLocTx;
    if isempty(ref_nodeLocTx)
        nodeLocTx{ii}.x=1.390;
        nodeLocTx{ii}.y=1.719;
        nodeLocTx{ii}.z=0.8;
        ref_nodeLocTx=[nodeLocTx{ii}.x,nodeLocTx{ii}.y,nodeLocTx{ii}.z];% reference TX just for scatterer deployment
    end
    for jj=1:config.num_RX
        mpc=[]; 
        channel_attempts=0;
        is_LOS=1;
        while((~config.Topology.channel_req || (config.Topology.channel_req && isempty(mpc))) && (UEgrid_idx <= Num_UE_grid_points))
            channel_attempts = channel_attempts + 1;   
            % Generate the UE location:
            nodeLocRx{jj}.x=UE_points(UEgrid_idx,1);
            nodeLocRx{jj}.y=UE_points(UEgrid_idx,2);
            if isfield(config.Topology,'UE_height')
                nodeLocRx{jj}.z = config.Topology.UE_height;
            else
                nodeLocRx{jj}.z = unifrnd(0.1,Dz);
            end
            UEgrid_idx=UEgrid_idx+1;
            % Compute the TX-RX distance
            % Run Raytracing Function and Generate Output
            % Input Node Related Parameters
            %nodeCfg = getNodeConfiguration(nodeLocTx{ii},nodeLocRx{jj});%exact nodes for ray tracing
            if isempty(ref_nodeLocRx)
                nodeLocRx{jj}.x=4.017;
                nodeLocRx{jj}.y=2.198;
                nodeLocRx{jj}.z=0.8;
                ref_nodeLocRx=[nodeLocRx{jj}.x,nodeLocRx{jj}.y,nodeLocRx{jj}.z];% reference RX just for scatterer deployment
            end
            ddist = sqrt((nodeLocTx{ii}.x-nodeLocRx{jj}.x)^2 +  (nodeLocTx{ii}.y-nodeLocRx{jj}.y)^2 +  (nodeLocTx{ii}.z-nodeLocRx{jj}.z)^2);
            nodeCfg = getNodeConfiguration(nodeLocTx{ii},nodeLocRx{jj},paraCfg);%exact nodes for ray tracing
            [PaaChannelInstances,scatterer_confg,LOS_op] = Raytracer(paraCfg, nodeCfg,ref_nodeLocTx,ref_nodeLocRx,scatterer_confg);
            if ~isempty(scatterer_confg)
                num=length(scatterer_confg(1,:));
            end
            for kk=1:paraCfg.numberOfTimeDivisions
                mpc = PaaChannelInstances{1,2,kk}; % node 1 = TX, node 2 = RX 
                if (isempty(mpc) || (size(mpc,1) == 1) && (mpc(1,9) < -120))
                    warning('No direct path or first-order reflection exists');
                    ampl_co1 = []; toa = []; az_tx = []; el_tx = []; az_rx = []; el_rx = []; 
                    mpc = [];phase =[];scatterer_confg(:,end-num+1:end)=[];ref_nodeLocRx=[];
                else
                    ampl_co1{kk}  = mpc(:,9); 
                    toa{kk}   = mpc(:,8); 
                    az_tx{kk} = mpc(:,10); 
                    el_tx{kk} = mpc(:,11);
                    az_rx{kk} = mpc(:,12); 
                    el_rx{kk} = mpc(:,13);
                    pl{kk}    = mpc(:,16);
                    ampl_co2{kk}   = mpc(:,17);
                    phase{kk} = mpc(:,18);
                    ampl_cx1{kk} = mpc(:,19);
                    doppler{kk} = mpc(:,20);
                    ampl_cx2{kk} = mpc(:,22);

                    [toa{kk}, id] = sort(toa{kk});

                    ampl_co1{kk}  = ampl_co1{kk}(id);
                    az_tx{kk} = az_tx{kk}(id); 
                    el_tx{kk} = el_tx{kk}(id); 
                    az_rx{kk} = az_rx{kk}(id); 
                    el_rx{kk} = el_rx{kk}(id);
                    pl{kk}    = pl{kk}(id);
                    ampl_co2{kk}   = ampl_co2{kk}(id);
                    phase{kk}= phase{kk}(id);
                    ampl_cx1{kk} = ampl_cx1{kk}(id);
                    doppler{kk}=doppler{kk}(id);
                    ampl_cx2{kk} = ampl_cx2{kk}(id);
                    
                    % Verify whether first path corresponds to LOS path by comparing time-of-arrival to toa of LOS distance
                    if LOS_op==1 && (config.Topology.channel_req~=3)
                        is_LOS = 1; disp('los path exist...');
                    elseif (config.Topology.channel_req==2) && LOS_op==0% If required LOS is missing, redo ray tracing.
                        is_LOS = 0;
                        mpc = [];
                        scatterer_confg(:,end-num+1:end)=[]; ref_nodeLocRx=[];disp('No los path exists');
                    elseif (config.Topology.channel_req==3) && LOS_op==1
                        is_LOS = 1;
                        mpc=[];
                        scatterer_confg(:,end-num+1:end)=[];ref_nodeLocRx=[];disp('los path exist...');
                    else 
                        is_LOS = 0;
                        disp('No los path exists');
                    end
                end
                if(config.Topology.channel_req == 0)
                    % It is ok to return empty channels, return
                    break;
                else
                    % Start over again?
                end
            end
        end
        s.Topology.UE_x(jj) = nodeLocRx{jj}.x;
        s.Topology.UE_y(jj) = nodeLocRx{jj}.y;
        s.Topology.UE_z(jj) = nodeLocRx{jj}.z;
        s.chan_config.Topology=s.Topology;
        % some intermediate plots:
        if config.is_plot
            figure(651);  
            subplot(2,1,1);  plot(nodeLocRx{jj}.x,nodeLocRx{jj}.y,'kX','MarkerSize',10,'MarkerFaceColor','k');
            xlabel('x (m)'); ylabel('y (m)'); 
            subplot(2,1,2);  plot3(nodeLocRx{jj}.x,nodeLocRx{jj}.y,nodeLocRx{jj}.z,'kX','MarkerSize',10,'MarkerFaceColor','k');
            xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); %grid minor;
        end
        % Normalize amplitudes: 
        % if isfield(config, 'useNormalizedPower') && (config.useNormalizedPower == 1)
        %     PL = 20 * log10(4 * pi * ddist * config.carrier_freq / 3e8);
        %     ampl = (10.^(-ampl/20)) / (10^(-PL/20)) ;  % Q++: perhaps we should rename this pathloss : this is magnitude squared (power)
        %     ampl = 20*log10(ampl);
        % end
        ampl = {ampl_co1,ampl_cx2;ampl_cx1,ampl_co2};
        output{ii,jj} = {ampl,toa, az_tx, el_tx, az_rx, el_rx,channel_attempts, ddist, is_LOS,pl,phase,doppler};
        % Q++: there are some warnings when I execute the code. can we check them together ?
    end                                                                                                                                                         
end
cd ../../..
end

%% GETSYSTEMPARAMETERS converts configuration parameters into system
% parameters used by the NIST ray tracer
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               %  Inputs:
%      config - configuration struct
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
%  Outputs:
%      param - struct with simulation paramete                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   rs
%
%  Author:  B. De Beelde <Brecht.DeBeelde@UGent.be>
%  Version: 1.2
%  September 2021; Last revision: June 17, 2022
function param = getSystemParameters(config)

    param = struct( ...
        'environmentFileName', 'conf_room.amf', ...
        'indoorSwitch', config.Topology.indoor, ...
        'numberOfTimeDivisions', config.Topology.timedivision, ...
        'referencePoint', [0 0 0], ...
        'selectPlanesByDist', inf, ...
        'switchDiffuseComponent', config.Topology.scattering, ...
        'switchQDModel', 'DBARC', ...
        'totalNumberOfReflections', config.Topology.reflectionorder, ...
        'totalTimeDuration', config.Topology.timeduration, ...
        'switchSaveVisualizerFiles', 1, ...
        'carrierFrequency', config.carrier_freq, ...
        'qdFilesFloatPrecision', 6, ...
        'useOptimizedOutputToFile', 1, ...
        'outputFormat', 'json', ...
        'materialLibraryPath', 'materialLibrary_conf.csv', ...
        'inputScenarioName', 'DBARC_CIR_generation', ...
        'diffusePathGainThreshold','-Inf',...
        'reflectionLoss', 10, ...
        'numberOfNodes', 2, ...
        'isInitialOrientationOn', 0, ...
        'isDeviceRotationOn', 0, ...
        'isPaaCentered', 1, ...
        'enablePhaseOutput',1,...
        'confg_dim', config.Topology.confg_dim, ...
        'scatterer_num', config.Topology.scatterer_number, ...
        'scatterer_dis',config.Topology.scatterer_distribution,...
        'scatterer_dim',config.Topology.scatterer_dimension,...
        'scatterer_height_max',config.Topology.scatterer_height_max,...
        'mat_confg',config.Topology.material_configuration,...
        'Ptx',config.Topology.transmit_power,...
        'antenna',config.Topology.antenna,...
        'polarization',config.Topology.polarization,...
        'scattering',config.Topology.scattering...
    );
end

%% GETNODECONFIGURATION creates node configurations used by the ray tracer, 
% based on a specific TX and RX antenna
%
%  Inputs:
%      nodeLocTx - struct with X,Y,Z coordinates of TX antenna
%      nodeLocRx - struct with X,Y,Z coordinates of RX antenna
%  Outputs:
%      nodeCfg - configuration struct used by NIST ray tracer
%
%  Author:  B. De Beelde <Brecht.DeBeelde@UGent.be>ad
%  Version: 1.2
%  September 2021; Last revision: June 17, 2022
function nodeCfg = getNodeConfiguration(nodeLocTx, nodeLocRx,paraCfg)

    nodeCfg = struct;
    nodeCfg.nodeLoc = ...
      [nodeLocTx.x nodeLocTx.y nodeLocTx.z;
       nodeLocRx.x nodeLocRx.y nodeLocRx.z];
    nodeCfg.nodeAntennaOrientation = cell(2, 1);
    nodeCfg.nodeAntennaOrientation{1,1} = ...
      [0 0 0];
    nodeCfg.nodeAntennaOrientation{2,1} = ...
      [0 0 0];
    nodeCfg.nodePosition = reshape([nodeLocTx.x nodeLocTx.y nodeLocTx.z nodeLocRx.x nodeLocRx.y nodeLocRx.z], paraCfg.numberOfTimeDivisions, 3, 2);
    nodeCfg.nodeRotation = reshape([0 0 0 0 0 0], 1, 3, 2);
    %tmp = zeros(paraCfg.numberOfTimeDivisions,3,2);
    %tmp(:,:,1)=[1,	19,	2.60000000000000;...
%3,	17,	2.60000000000000;...
%5,	15,	2.60000000000000;....
% 7,	13,	2.60000000000000;...
% 9,	11,	2.60000000000000;...
% 11,	9,	2.60000000000000;...
% 13,	7,	2.60000000000000;...
% 15,	5,	2.60000000000000;...
% 17,	3,	2.60000000000000;...
% 19,	1,	2.60000000000000];
%     tmp(:,:,2)=[2.26358981136190,	11.1822143065056,	1.60000000000000;...
% 2.61939541780178,	11.5656815514044,	1.60000000000000;...
% 3.02837979677732,	11.8918357101060,	1.60000000000000;...
% 3.48140691308855,	12.1533910376585,	1.60000000000000;...
% 3.96835688956183,	12.3445048202515,	1.60000000000000;...
% 4.47835206838339,	12.4609078917880,	1.60000000000000;...
% 5.52164793161661,	12.4609078917880,	1.60000000000000;...
% 6.03164311043816,	12.3445048202515,	1.60000000000000;...
% 6.51859308691145,	12.1533910376585,	1.60000000000000;...
% 6.97162020322268,	11.8918357101060,	1.60000000000000];
    %nodeCfg.nodePosition = tmp;
    nodeCfg.nodeRotation = zeros(paraCfg.numberOfTimeDivisions,3,2);
    nodeCfg.paaInfo = cell(1, 2);
    nodeCfg.paaInfo{1} = struct;
    nodeCfg.paaInfo{1}.nPaa = 1;
    nodeCfg.paaInfo{1}.centroids = 1;
    nodeCfg.paaInfo{1}.paaInCluster = cell(1, 1);
    nodeCfg.paaInfo{1}.paaInCluster{1} = 1;
    nodeCfg.paaInfo{1}.centroidsShift = cell(1, 1);
    nodeCfg.paaInfo{1}.centroidsShift{1} = ...
      [0 0 0];
    nodeCfg.paaInfo{1}.PAA_loc = reshape( ...
      [nodeLocTx.x nodeLocTx.y nodeLocTx.z], 1, 3, 1);
%     tmp=zeros(paraCfg.numberOfTimeDivisions,1,3);
%     tmp(:,1,:)=[1,	19,	2.60000000000000;...
% 3,	17,	2.60000000000000;...
% 5,	15,	2.60000000000000;....
% 7,	13,	2.60000000000000;...
% 9,	11,	2.60000000000000;...
% 11,	9,	2.60000000000000;...
% 13,	7,	2.60000000000000;...
% 15,	5,	2.60000000000000;...
% 17,	3,	2.60000000000000;...
% 19,	1,	2.60000000000000];
    %nodeCfg.paaInfo{1}.PAA_loc=tmp;
    nodeCfg.paaInfo{1}.orientation = cell(1, 1);
    nodeCfg.paaInfo{1}.orientation{1} = ...
      [0 0 0];
    nodeCfg.paaInfo{1}.generationMethod = 2;
    nodeCfg.paaInfo{1}.nodePAAInfo = cell(1, 1);
    nodeCfg.paaInfo{1}.nodePAAInfo{1} = struct;
    nodeCfg.paaInfo{1}.nodePAAInfo{1}.centroid_id = 1;
    nodeCfg.paaInfo{1}.nodePAAInfo{1}.tot_channel = 1;
    nodeCfg.paaInfo{1}.nodePAAInfo{1}.indep_stoch_channel = 1;
    nodeCfg.paaInfo{1}.nodePAAInfo{1}.rotated_channel = 1;
    nodeCfg.paaInfo{1}.nodePAAInfo{1}.paa_id = 1;
    nodeCfg.paaInfo{1}.centroidTimePosition = reshape( ...
      [nodeLocTx.x nodeLocTx.y nodeLocTx.z], paraCfg.numberOfTimeDivisions, 3, 1);
%     tmp=zeros(paraCfg.numberOfTimeDivisions,1,3);
%     tmp(:,1,:)=[1,	19,	2.60000000000000;...
% 3,	17,	2.60000000000000;...
% 5,	15,	2.60000000000000;....
% 7,	13,	2.60000000000000;...
% 9,	11,	2.60000000000000;...
% 11,	9,	2.60000000000000;...
% 13,	7,	2.60000000000000;...
% 15,	5,	2.60000000000000;...
% 17,	3,	2.60000000000000;...
% 19,	1,	2.60000000000000];
    %nodeCfg.paaInfo{1}.centroidTimePosition=tmp;
    nodeCfg.paaInfo{1}.nPAA_centroids = 1;
    nodeCfg.paaInfo{1}.node_centroid = [nodeLocTx.x nodeLocTx.y nodeLocTx.z];
%     nodeCfg.paaInfo{1}.node_centroid = [1,	19,	2.60000000000000;...
% 3,	17,	2.60000000000000;...
% 5,	15,	2.60000000000000;....
% 7,	13,	2.60000000000000;...
% 9,	11,	2.60000000000000;...
% 11,	9,	2.60000000000000;...
% 13,	7,	2.60000000000000;...
% 15,	5,	2.60000000000000;...
% 17,	3,	2.60000000000000;...
% 19,	1,	2.60000000000000];
    nodeCfg.paaInfo{2} = struct;
    nodeCfg.paaInfo{2}.nPaa = 1;
    nodeCfg.paaInfo{2}.centroids = 1;
    nodeCfg.paaInfo{2}.paaInCluster = cell(1, 1);
    nodeCfg.paaInfo{2}.paaInCluster{1} = 1;
    nodeCfg.paaInfo{2}.centroidsShift = cell(1, 1);
    nodeCfg.paaInfo{2}.centroidsShift{1} = ...
      [0 0 0];
    nodeCfg.paaInfo{2}.PAA_loc = reshape( ...
      [nodeLocRx.x nodeLocRx.y nodeLocRx.z], 1, 3, 1);
%     tmp(:,1,:)=[18,	19,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 10,	17,	1.60000000000000;...
% 10,	10,	1.60000000000000;...
% 5,	10,	1.60000000000000;...
% 3,	1,	1.60000000000000;...
% 1,	1,	1.60000000000000];

    %nodeCfg.paaInfo{2}.PAA_loc=tmp;
    nodeCfg.paaInfo{2}.orientation = cell(1, 1);
    nodeCfg.paaInfo{2}.orientation{1} = ...
      [0 0 0];
    nodeCfg.paaInfo{2}.generationMethod = 2;
    nodeCfg.paaInfo{2}.nodePAAInfo = cell(1, 1);
    nodeCfg.paaInfo{2}.nodePAAInfo{1} = struct;
    nodeCfg.paaInfo{2}.nodePAAInfo{1}.centroid_id = 1;
    nodeCfg.paaInfo{2}.nodePAAInfo{1}.tot_channel = 1;
    nodeCfg.paaInfo{2}.nodePAAInfo{1}.indep_stoch_channel = 1;
    nodeCfg.paaInfo{2}.nodePAAInfo{1}.rotated_channel = 1;
    nodeCfg.paaInfo{2}.nodePAAInfo{1}.paa_id = 1;
    nodeCfg.paaInfo{2}.centroidTimePosition = reshape( ...
      [nodeLocRx.x nodeLocRx.y nodeLocRx.z],paraCfg.numberOfTimeDivisions,3, 1);
    %tmp=zeros(paraCfg.numberOfTimeDivisions,1,3);
%     tmp(:,1,:)=[18,	19,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 10,	17,	1.60000000000000;...
% 10,	10,	1.60000000000000;...
% 5,	10,	1.60000000000000;...
% 3,	1,	1.60000000000000;...
% 1,	1,	1.60000000000000];

    %nodeCfg.paaInfo{2}.centroidTimePosition=tmp;
    nodeCfg.paaInfo{2}.nPAA_centroids = 1;
    nodeCfg.paaInfo{2}.node_centroid = [nodeLocRx.x nodeLocRx.y nodeLocRx.z];
%     nodeCfg.paaInfo{2}.node_centroid= [18,	19,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 18,	18,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 15,	14,	1.60000000000000;...
% 10,	17,	1.60000000000000;...
% 10,	10,	1.60000000000000;...
% 5,	10,	1.60000000000000;...
% 3,	1,	1.60000000000000;...
% 1,	1,	1.60000000000000];

end