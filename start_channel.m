%
% This script launches a simulation of high-throughput sub-Thz communications for advanced radios.
%
% This software is copyrighted by Imec.
%
% Some guideline on the tool configuration:
%===========================================
% 1: To investigate the influence on scatterer distribution on the channel
%    modelling, it's better to set the 'scatterer_dimension' to 0.
%
% 2: When configuring large scatterers, the 'scatterer_number' shouldn't
% be high, better less than 10.
%
% 3: The 'confg_dim' should to be larger than [5 5 2], otherwise the
% scatterer deployment will be affected.
%
% 4: 'AP_arch' is better to be set to 'ceil-random', and 'UE_arch' is
% better to be set to 'random'.
%
% 5: exact scatterer dimension: small (0.25*0.25*0.125)~(0.5*0.5*0.25),
% medium (1*1*0.5)~(2*2*1), large (2.5*2.5*1.25)~(5*5*2.5), mixed
% (combination of the above three dimensions)
%
% 6: When getting the error "inapproiate configuration", 'confg_dim',
% 'scatterer_dim',or 'scatterer_num' should be adjusted.
%
% 7: The scatterer distribution configuration should be an array.
% --------------------------------------------------------------
% There three different configurations:
% **7.1) Disk (2D, including ring and circle): [0, distribution, Radius_inner, Radius_outer, mean, variance]
%       - 'distribution' of the scatterer location:
%                0 = uniform distribution
%                1 = hyperbolic distribution
%                2 = Gaussian distribution
%                3 = Rayleigh distribution
%       - 'Radius_inner' and  'Radius_outer' are resp. inner and outer circle radius
%                Ring   if Radius_inner > 0
%                Circle if Radius_inner = 0
%           Notes for Statistical distribution
%       - 'mean' and 'variance' of the distributions. The mean must be within the inner radius and outer radius.
%               Note for Rayleigh and Hyperbolic distribution,'mean' is not used and only 'variance' is used.
%               For uniform distribution, the mean and var are not used.
%
% **7.2) Sphere (3D): [1, distribution, Radius_inner, Radius_outer, mean, variance]
%       - 'distribution' of the scatterer location:
%                0 = uniform distribution
%                1 = hyperbolic distribution
%                2 = Gaussian distribution
%                3 = Rayleigh distribution
%       - 'Radius_inner' and  'Radius_outer' are resp. inner and outer sphere radius
%       Notes for Statistical distribution
%       - 'mean' and 'variance' of the distributions. The mean must be within the inner radius and outer radius.
%               Note for Rayleigh and Hyperbolic distribution,'mean' is not used and only 'variance' is used.
%               For uniform distribution, the mean and var are not used.
% **7.3) Cluster: [2, Number_Of_Clusters,Max_Number_Of_Scatter_per_Clusters,Cluster_Maximum_Radius,...
%                                      cluster_center_distribution, scatterer_in_cluster_distribution]
%      - 'cluster_center_distribution': the statistical distribution of cluster centers:
%                0 = uniform distribution
%                1 = hyperbolic distribution (std=1)
%                2 = Gaussian distribution (mean=center of the dimension
%                (0.5*confg_dim(1),0.5*confg_dim(2),0.5*confg_dim(3), std=5)
%                3 = Rayleigh distribution (std=1)
%      - 'scatterer_in_cluster_distribution': the statistical distribution of scatterers in a cluster
%         (The geometrical distribution of scatterers in a cluster is assumed to be a sphere)
%                0 = uniform distribution
%                1 = hyperbolic distribution (std=1)
%                2 = Gaussian distribution (mean=center of the dimension
%                (0.5*confg_dim(1),0.5*confg_dim(2),0.5*confg_dim(3), std=5)
%                3 = Rayleigh distribution (std=1)
%       Cluster details can be found in the raytracer/mycluster.m file
% **7.4) Uniform (Geometrical): the whole valid space will be deployed with scatterers, which
%       selects certain configured number of scatterer grid porints
%       randomly.

%  ===> Q++: we wil need also to be able to pass some static configurations to the Dband engine
%  ===> Q++: return the scatterer grid points or at least their number in case of uniform scatterers

%
% Further remarks:
% ----------------
%
% 8.1) **For disk and sphere distribution, the scatterers are deployed around
%   a certain reference UE point (in x,y plane). The height of the
%   scatterers are determined by min(UE_height, scatterer_height_max)
%      **For cluster distribution, the scatterers are deployed in the whole confg_dimension
% 8.2) The scatterer configurations are returned in scatterer_confg in which
%       - the first three rows are the scatterer center coordinates (x,y,z), and
%       - the last three rows are the scatterer length, width, height (x,y,z).
%       - the number of columns are the scatterer number
%
%
% 9) distributed MIMO can be reliazed by adjusting the Kue and Kap, but
% current AP and UE selection algorithm just supports single AP.
%
% 10) indoor/outdoor scenarios can also be configured.
%
%
%
%
% Modified by Yunqi Feng<Yunqi.Feng@UGent.be>
% December 2023
%clear all ; close all; clc

% addpath classes demo utilities data
addpath(genpath('functions')) % functional call also including sub-folders

global s
global num_RX
global num_TX
sample_rate = 10; % GHz
theta = 0:1:360; % degree
theta_1 = -90:1:90; % degree

tau = 0:0.1:200;    % ns

%% Parameter initialization
s = struct(... % scenario parameters
    ...%%% Channel and simulation length control
    'random_seed',      1,    ...   % initial random seed value (if >= 0 only, otherwise not initialized)
    'num_channels',     1,   ...   % number of channel instances tested
    ...% channel model
    'carrier_freq',     140e9,  ... % carrier frequency, Hz
    'BW',               10e9, ...   % RF bandwidth, Hz (symmetrical around the carrier)
    'Topology', struct( ...
    'indoor',         1, ...% 1=indoor, 0=outdoor
    'reflectionorder',  2,...% maximum reflection order considered in the simulation, either 1 or 2
    'channel_req',    1, ... % 0 = no constraints, 1 = channel should be non-empty, 2 = LOS path required
    'AP_arch',        'ceil-collocated', ...  % AP architecture: 'ceil-piazza' , 'ceil-collocated', 'ceil-grid' , 'ceil-random': just a placeholer for future developpements
    'AP_gridNum',     [5 5],   ... % number of points in the AP grid along x/y
    'AP_safty',       [0.5 0.5], ... % the safety distance wrt the walls when placing the AP grid points
    'AP_height',      25+25*rand,   ...      % The height of the AP antennas
    'UE_arch',        'random', ...  % UE architecture: 'grid' , 'random'
    'UE_gridNum',     [2 2],   ... % number of points in the UE grid along x/y
    'UE_safty',       [1.5 1.5], ... % the safety distance wrt the walls when placing the UE grid points
    'UE_height',      1.5+21*rand,   ...      % The height of the AP antennas e.g. 1.65 , 0.5
    'confg_dim',      [100 100 100],  ...  % (x,y,z): dimensions of the room (rectangular, meters)
    'transmit_power', 0,... % transmit power in dBm, in the modeling, the noise level is -120 dBm
    'antenna',        'scan',... % TX antenna configuration-RX antenna configuration: 'omni'/'dir'/'scan' (dir: directional antenna pattern;scan: a channel sounding approach that obtains omni-directional channel characteristics, high gain. to model omnidirectional channel, it's better to choose 'scan' option than 'omni')
    'polarization',   'V-V',...,%TX antenna polarization-RX antenna polarization: 'V-V','H-H','dual'
    'scattering',   0,...,% 0 = without scattering, 1 = with scattering. suggestion: set to 0 because currently only first-order scattering can be modeled and the contribution is quite small if the material surfaces are not very rough
    'timeduration',   1, .... % total time evolution (used for doppler modelling, currently not used)
    'timedivision',   1, ... % the number of time samples to simulate(used for doppler modelling, currently not used)
    'scatterer_number',       16, ...% the total number of scatterers 
    'scatterer_dimension',    2,...%the dimension of the scatterer,0=small,1=medium,2=large,3=mixed
    'scatterer_height_max',   100,...% maximum scatterer height (meter) (cannot exceed AP height)
    'scatterer_distribution', 3, ...% See section 7 in 'dependent_parameters' e.g. [0, distribution, Radius_inner, Radius_outer, mean, variance]
    'material_configuration',[randi([1 2]) randi([1 3]) randi([1 2]) randi([1 3]) 1 4] ...% for indoor:material of front/back/left/right wall/ceiling/floor. 1=plaster,2=brick, 3=glass, 4=wood, and ceiling and floor can not be configured by glass. Scatterer materials are configured inside the ray tracer
    ) ...
    );
%  ===> Q++: It would be good to add some sanity check son the scatterer configurations
% +20*rand 9 randi([3 9])
% outdoor [80+50*rand 80+50*rand 50] randi([6,9]) 2.5+20*rand 1.5+5*rand
% indoor [40+10*rand 40+10*rand 5] randi([6,9]) 2.5 1.5
% sanity checks:
% if s.Topology.reflectionorder~=1 && s.Topology.reflectionorder~=2
%     error('maximum reflection order considered in the simulation can just be 1 or 2');
% end
if (s.Topology.AP_height > s.Topology.confg_dim(3)) && s.Topology.indoor==1
    disp('AP antennas height must be smaller that the  room height ==>  will be adjusted accordingly');
    s.Topology.AP_height = s.Topology.confg_dim(3);
end
if s.Topology.UE_height > s.Topology.confg_dim(3)&&s.Topology.indoor==1
    disp('UE antennas height must be smaller than room height ==>  will be adjusted accordingly');
    s.Topology.UE_height = s.Topology.confg_dim(3);
end
if s.Topology.UE_height > s.Topology.AP_height
    disp('UE antennas height must be smaller than AP height ==>  will be adjusted accordingly');
    s.Topology.UE_height = s.Topology.AP_height;
end
%if s.Topology.scatterer_height_max > s.Topology.AP_height
    %disp('scatterer height must be smaller than AP height ==>  will be adjusted accordingly');
    %s.Topology.scatterer_height_max = s.Topology.UE_height+rand*(s.Topology.AP_height-s.Topology.UE_height);
%end
if length(s.Topology.material_configuration)~=6 && s.Topology.indoor==1
    error('Four materials should be chosen for indoor scenarios')
elseif s.Topology.material_configuration(5)==3 || s.Topology.material_configuration(6)==3
    error('ceiling and floor can not be configured by glass')
end
if ~isempty(s.Topology.material_configuration) && s.Topology.indoor==0
    %disp('Outdoor material characerization have been configured inside.')
end

% Define the numbers of TXs and RXs
if ~isfield(s,'Kue')
    s.Kue = 1;
end
if ~isfield(s,'Kap')
    s.Kap=1;
end
% Generate the grid points for the APs and UEs
generate_APgrid(s.Kap);
generate_UEgrid(s.Kue);
% Some assignments:
s.Dband_chConfig = struct('Topology', s.Topology, ...
    'carrier_freq', s.carrier_freq, 'num_TX', 1, 'num_RX', s.Kue, ...
    'gain', struct('tx', 1, 'rx', 1), 'txpower', 1, 'useNormalizedPower', 0, ... % For for the time being: 'gain'/'power' are not relevant/used while 'useNormalizedPower' must be set to 0
    'is_plot',0);
if isfield(s,'paraCfg')
    s = rmfield(s,'paraCfg');
end

T=s.Dband_chConfig.num_TX;
R=s.Dband_chConfig.num_RX;
%rng(33987);
rng('shuffle')
%% channel construction:
%
for ch = 1:1 %s.num_channels
    
    if isfield(s.Dband_chConfig.Topology,'nodeLocTx'), s.Dband_chConfig.Topology = rmfield(s.Dband_chConfig.Topology,'nodeLocTx'); end
    if isfield(s.Dband_chConfig.Topology,'nodeLocRx'), s.Dband_chConfig.Topology = rmfield(s.Dband_chConfig.Topology,'nodeLocRx'); end
    
    % Generate the geometrical channel
    [output,scatterer_confg]  = Dband_model_engine();
    
    % Parse the data:
    s.Dband_scatterer = scatterer_confg; clear scatterer_confg
    clear ampl phas
    for ii=1:s.Dband_chConfig.num_TX
        for jj=1:s.Dband_chConfig.num_RX
            ampl                           = output{ii,jj}{1};  %  (dB): gains  of the propagation paths
            s.Dband_channel{ii,jj}.ampl_co1  =ampl{1,1};
            ampl_co1                        = s.Dband_channel{ii,jj}.ampl_co1{1,1};
            s.Dband_channel{ii,jj}.ampl_co2  =ampl{2,2};
            ampl_co2                        = s.Dband_channel{ii,jj}.ampl_co2{1,1};
            s.Dband_channel{ii,jj}.ampl_cx1  =ampl{2,1};
            ampl_cx1                        = s.Dband_channel{ii,jj}.ampl_cx1{1,1};
            s.Dband_channel{ii,jj}.ampl_cx2  =ampl{1,2};
            ampl_cx2                        = s.Dband_channel{ii,jj}.ampl_cx2{1,1};
            s.Dband_channel{ii,jj}.toa     = output{ii,jj}{2}{1,1};  % (ns): time of arrival
            s.Dband_channel{ii,jj}.az_tx   = output{ii,jj}{3}{1,1};
            s.Dband_channel{ii,jj}.el_tx   = output{ii,jj}{4}{1,1};  % ( ): AoD elevation
            s.Dband_channel{ii,jj}.az_rx   = output{ii,jj}{5}{1,1};  % ( ): AoA azimut
            s.Dband_channel{ii,jj}.el_rx   = output{ii,jj}{6}{1,1};  % ( ): AoA elevation
            s.Dband_channel{ii,jj}.attempts= output{ii,jj}{7}; % Number of attempts for a successful covered... can be used for black hole coverage spots analysis
            s.Dband_channel{ii,jj}.ddist   = output{ii,jj}{8};% Tx-RX distance
            s.Dband_channel{ii,jj}.is_LOS  = output{ii,jj}{9}; % whether or not a LOS path is present
            phas                           = output{ii,jj}{11}{1,1}; % (rad): phases  of the propagation paths (see also output{ii,jj}{1} for the gains)
            s.Dband_channel{ii,jj}.pl      = output{ii,jj}{10}{1,1};
            s.Dband_channel{ii,jj}.phas  = output{ii,jj}{11}{1,1}; % (rad): phases  of the propagation paths (see also output{ii,jj}{1} for the gains)
            s.Dband_channel{ii,jj}.doppler  = output{ii,jj}{12}{1,1}; % (rad): phases  of the propagation paths (see also output{ii,jj}{1} for the gains)
            if strcmp(s.Topology.polarization, 'dual')
                s.Dband_channel{ii,jj}.cmplxGain_co1  = 10.^(ampl_co1/20) .* exp(1i*2*pi*phas);
                s.Dband_channel{ii,jj}.cmplxGain_co2  = 10.^(ampl_co2/20) .* exp(1i*2*pi*phas);
                s.Dband_channel{ii,jj}.cmplxGain_cx1  = 10.^(ampl_cx1/20) .* exp(1i*2*pi*phas);
                s.Dband_channel{ii,jj}.cmplxGain_cx2  = 10.^(ampl_cx2/20) .* exp(1i*2*pi*phas);
            elseif strcmp(s.Topology.polarization, 'V-V') || strcmp(s.Topology.polarization, 'H-H')
                s.Dband_channel{ii,jj}.cmplxGain_co1  = 10.^(ampl_co1/20) .* exp(1i*2*pi*phas);
                s.Dband_channel{ii,jj}.cmplxGain_cx1  = 10.^(ampl_cx1/20) .* exp(1i*2*pi*phas);
            end
            if strcmp(s.Topology.antenna,'scan')
                if s.Dband_channel{ii,jj}.is_LOS == 0
                    az_tx = s.Dband_channel{ii,jj}.az_tx;
                    el_tx = s.Dband_channel{ii,jj}.el_tx;
                    az_rx = s.Dband_channel{ii,jj}.az_rx;
                    el_rx = s.Dband_channel{ii,jj}.el_rx;
                    ampl_co1 = s.Dband_channel{ii,jj}.ampl_co1{1,1};
                    sum_ampl = [];
                    for k = 1:length(s.Dband_channel{1,1}.az_tx)
                        ampl_co1_adjusted = ampl_co1(k);
                        ini_az_tx = az_tx(k);
                        ini_el_tx = el_tx(k);
                        ini_ampl_co1 = ampl_co1(k);
                        % Find the index of the first occurrence of the extracted value
                        idx_az_tx = find(az_tx == ini_az_tx, 1);
                        idx_el_tx = find(el_tx == ini_el_tx, 1);
                        idx_ampl_co1 = find(ampl_co1 == ini_ampl_co1, 1);
                        
                        % Remove only the first occurrence of the extracted value
                        rest_az_tx = az_tx;
                        rest_az_tx(idx_az_tx) = [];  % Remove the first occurrence
                        
                        rest_el_tx = el_tx;
                        rest_el_tx(idx_el_tx) = [];  % Remove the first occurrence
                        
                        rest_ampl_co1 = ampl_co1;
                        rest_ampl_co1(idx_ampl_co1) = [];  % Remove the first occurrence

                        rel_az_tx = rest_az_tx - ini_az_tx;
                        rel_el_tx = rest_el_tx - ini_el_tx;
                        ini_az_rx = az_rx(k);
                        ini_el_rx = el_rx(k);
                        ini_ampl_co1 = ampl_co1(k);
                        % Find the index of the first occurrence of the extracted value
                        idx_az_rx = find(az_rx == ini_az_rx, 1);
                        idx_el_rx = find(el_rx == ini_el_rx, 1);
                        idx_ampl_co1 = find(ampl_co1 == ini_ampl_co1, 1);
                        
                        % Remove only the first occurrence of the extracted value
                        rest_az_rx = az_rx;
                        rest_az_rx(idx_az_rx) = [];  % Remove the first occurrence
                        
                        rest_el_rx = el_rx;
                        rest_el_rx(idx_el_rx) = [];  % Remove the first occurrence
                        
                        rest_ampl_co1 = ampl_co1;
                        rest_ampl_co1(idx_ampl_co1) = [];  % Remove the first occurrence

                        rel_az_rx = rest_az_rx - ini_az_rx;
                        rel_el_rx = rest_el_rx - ini_el_rx;
                        for kk = 1:length(rest_az_tx)
                            G_t = horn_gain(rel_el_tx(kk),rel_az_tx(kk)) - 23.7;
                            G_r = horn_gain(rel_el_rx(kk),rel_az_rx(kk)) - 23.7;
                            ampl_co1_adjusted_tmp = rest_ampl_co1(kk) + G_t + G_r;
                            if ampl_co1_adjusted_tmp<-120
                                ampl_co1_adjusted_tmp = [];
                            else
                                ampl_co1_adjusted = [ampl_co1_adjusted ampl_co1_adjusted_tmp];
                            end
                        end
                        sum_ampl = [sum_ampl sum(10.^(ampl_co1_adjusted./10))];
                    end
                    [max_val, idx] = max(sum_ampl);
                    best_nlos_dir = 10 * log10(sum_ampl(idx));
                elseif s.Dband_channel{ii,jj}.is_LOS == 1
                    az_tx = s.Dband_channel{ii,jj}.az_tx;
                    el_tx = s.Dband_channel{ii,jj}.el_tx;
                    az_rx = s.Dband_channel{ii,jj}.az_rx;
                    el_rx = s.Dband_channel{ii,jj}.el_rx;
                    ampl_co1 = s.Dband_channel{ii,jj}.ampl_co1{1,1};
                    az_tx = az_tx(2:end);
                    el_tx = el_tx(2:end);
                    az_rx = az_rx(2:end);
                    el_rx = el_rx(2:end);
                    sum_ampl = [];
                    for k = 1:length(az_tx)
                        ampl_co1_adjusted = ampl_co1(k);
                        ini_az_tx = az_tx(k);
                        ini_el_tx = el_tx(k);
                        ini_ampl_co1 = ampl_co1(k);
                        % Find the index of the first occurrence of the extracted value
                        idx_az_tx = find(az_tx == ini_az_tx, 1);
                        idx_el_tx = find(el_tx == ini_el_tx, 1);
                        idx_ampl_co1 = find(ampl_co1 == ini_ampl_co1, 1);
                        
                        % Remove only the first occurrence of the extracted value
                        rest_az_tx = az_tx;
                        rest_az_tx(idx_az_tx) = [];  % Remove the first occurrence
                        
                        rest_el_tx = el_tx;
                        rest_el_tx(idx_el_tx) = [];  % Remove the first occurrence
                        
                        rest_ampl_co1 = ampl_co1;
                        rest_ampl_co1(idx_ampl_co1) = [];  % Remove the first occurrence

                        rel_az_tx = rest_az_tx - ini_az_tx;
                        rel_el_tx = rest_el_tx - ini_el_tx;
                        ini_az_rx = az_rx(k);
                        ini_el_rx = el_rx(k);
                        ini_ampl_co1 = ampl_co1(k);
                        % Find the index of the first occurrence of the extracted value
                        idx_az_rx = find(az_rx == ini_az_rx, 1);
                        idx_el_rx = find(el_rx == ini_el_rx, 1);
                        idx_ampl_co1 = find(ampl_co1 == ini_ampl_co1, 1);
                        
                        % Remove only the first occurrence of the extracted value
                        rest_az_rx = az_rx;
                        rest_az_rx(idx_az_rx) = [];  % Remove the first occurrence
                        
                        rest_el_rx = el_rx;
                        rest_el_rx(idx_el_rx) = [];  % Remove the first occurrence
                        
                        rest_ampl_co1 = ampl_co1;
                        rest_ampl_co1(idx_ampl_co1) = [];  % Remove the first occurrence

                        rel_az_rx = rest_az_rx - ini_az_rx;
                        rel_el_rx = rest_el_rx - ini_el_rx;
                        for kk = 1:length(rest_az_tx)
                            G_t = horn_gain(rel_el_tx(kk),rel_az_tx(kk)) - 23.7;
                            G_r = horn_gain(rel_el_rx(kk),rel_az_rx(kk)) - 23.7;
                            ampl_co1_adjusted_tmp = rest_ampl_co1(kk) + G_t + G_r;
                            if ampl_co1_adjusted_tmp<-120
                                ampl_co1_adjusted_tmp = [];
                            else
                                ampl_co1_adjusted = [ampl_co1_adjusted ampl_co1_adjusted_tmp];
                            end
                        end
                        sum_ampl = [sum_ampl sum(10.^(ampl_co1_adjusted./10))];
                    end
                    [max_val, idx] = max(sum_ampl);
                    best_nlos_dir = 10 * log10(sum_ampl(idx));
                end
            else
                best_nlos_dir = NaN;
                            
            end


        end
    end
    clear ii jj
    % Investigate parameters from transceivers of interest
    %Dband_plots();%only one co-polarization components are plotted but straightforward to plot other polarization components
end

%% co-polarized channel statistics
ampl_co_dB = ampl{1,1}{1,1};
ampl_cx_dB = s.Dband_channel{1,1}.ampl_cx1{1,1};
ampl_cx_dB = ampl_cx_dB(~isnan(ampl_cx_dB));
toa=s.Dband_channel{1,1}.toa ;
az_rx=s.Dband_channel{1,1}.az_rx;
el_rx=s.Dband_channel{1,1}.el_rx;
az_tx=s.Dband_channel{1,1}.az_tx;
el_tx=s.Dband_channel{1,1}.el_tx;
ddist=s.Dband_channel{1,1}.ddist;
ampl_co=10.^(ampl_co_dB/10);
ampl_cx=10.^(ampl_cx_dB/10);
pl_co_dB=s.Dband_channel{1,1}.pl;
pl_co=10.^(pl_co_dB/10);
num_ray=length(ampl_co);
toa=toa*1e9;
if s.Topology.channel_req==2
    %if length(pl_co)>1
        % pl_co_los=pl_co(1);
        % pl_co_nlos=sum(pl_co)-pl_co_los;
    pl_co_los = sum(pl_co);
    pl_co_nlos=NaN;
elseif s.Topology.channel_req==3
    pl_co_nlos=sum(pl_co);
    pl_co_los=NaN;
    %pl_co_los=pl_co(1);
elseif s.Topology.channel_req==1
    pl_co_los=NaN;
    pl_co_nlos=NaN;
end
toa_mean=sum(toa.*ampl_co)/sum(ampl_co);
toa_ds=sqrt(sum(toa.^2.*ampl_co)/sum(ampl_co)-toa_mean^2);
az_tx_mean=sum(exp(1*j*(az_tx*pi/180)).*ampl_co)/sum(ampl_co);
az_tx_as=180/pi*sqrt(sum(((cosd(az_tx)-real(az_tx_mean)).^2+(sind(az_tx)-imag(az_tx_mean)).^2).*ampl_co)/sum(ampl_co));
az_rx_mean=sum(exp(1*j*(az_rx*pi/180)).*ampl_co)/sum(ampl_co);
az_rx_as=180/pi*sqrt(sum(((cosd(az_rx)-real(az_rx_mean)).^2+(sind(az_rx)-imag(az_rx_mean)).^2).*ampl_co)/sum(ampl_co));
el_tx_mean=sum(exp(1*j*(el_tx*pi/180)).*ampl_co)/sum(ampl_co);
el_tx_as=180/pi*sqrt(sum(((cosd(el_tx)-real(el_tx_mean)).^2+(sind(el_tx)-imag(el_tx_mean)).^2).*ampl_co)/sum(ampl_co));
el_rx_mean=sum(exp(1*j*(el_rx*pi/180)).*ampl_co)/sum(ampl_co);
el_rx_as=180/pi*sqrt(sum(((cosd(el_rx)-real(el_rx_mean)).^2+(sind(el_rx)-imag(el_rx_mean)).^2).*ampl_co)/sum(ampl_co));

%clear output ampl_co az_tx el_tx az_rx el_rx attempts is_LOS
%clear output attempts is_LOS
function G = horn_gain(elevation, azimuth)
% Load the gain matrix from the .mat file
% currentDir = fileparts(mfilename('fullpath'));
% parentDir = fileparts(currentDir); % Get the parent directory
% Construct the full path to the data file
%dataFile = fullfile(parentDir, 'horn_gain.mat'); % Replace with the actual filename
dataFile = 'C:\Users\yufeng\OneDrive - UGent\WAVES-D-band-channel-models\WAVES-D-band-Channel-main_Nov\functions\Channel_model\DBARC\\horn_gain.mat';
gain_matrix = load(dataFile).G;
[n_azimuth, n_elevation] = size(gain_matrix);

% Define the grid for elevation and azimuth angles
%elevation_angles = linspace(-90, 90, n_elevation); % Elevation from -90 to 90 degrees
%azimuth_angles = linspace(-180, 180, n_azimuth);      % Azimuth from 0 to 360 degrees

elevation_angles = linspace(0, 180, n_elevation); % Elevation from -90 to 90 degrees
azimuth_angles = linspace(0, 360, n_azimuth);      % Azimuth from 0 to 360 degrees

% Find the closest indices for the input angles
[~, elevation_idx] = min(abs(elevation_angles - elevation));
[~, azimuth_idx] = min(abs(azimuth_angles - azimuth));

% Select the antenna gain
G = gain_matrix(azimuth_idx, elevation_idx);
end



