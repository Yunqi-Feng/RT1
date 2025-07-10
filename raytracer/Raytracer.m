function [outputPaaTime,scatterer_confg,isLos] = Raytracer(paraCfgInput, nodeCfgInput,ref_nodeLocTx,ref_nodeLocRx,scatterer_confg)
%%RAYTRACER generates the QD channel model.
% Inputs:
% paraCfgInput - Simulation configuration
% nodeCfgInput - Node configuration 

% Outputs:
% outputPaaTime: PAA channel output
% scatterer_confg: configured scatterers in the floor map described by
% their centers and dimensions

%--------------------------Software Disclaimer-----------------------------
%
% NIST-developed software is provided by NIST as a public service. You may 
% use, copy and distribute copies of the software in any medium, provided 
% that you keep intact this entire notice. You may improve, modify and  
% create derivative works of the software or any portion of the software, 
% and you  may copy and distribute such modifications or works. Modified 
% works should carry a notice stating that you changed the software and  
% should note the date and nature of any such change. Please explicitly  
% acknowledge the National Institute of Standards and Technology as the 
% source of the software.
% 
% NIST-developed software is expressly provided "AS IS." NIST MAKES NO
% WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY OPERATION  
% OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTY OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT AND 
% DATA ACCURACY. NIST NEITHER REPRESENTS NOR WARRANTS THAT THE OPERATION OF 
% THE SOFTWARE WILL BE UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS 
% WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS  
% REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF, INCLUDING BUT 
% NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR USEFULNESS OF 
% THE SOFTWARE.
%
% You are solely responsible for determining the appropriateness of using
% and distributing the software and you assume all risks associated with  
% its use, including but not limited to the risks and costs of program 
% errors, compliance with applicable laws, damage to or loss of data, 
% programs or equipment, and the unavailability or interruption of 
% operation. This software is not intended to be used in any situation  
% where a failure could cause risk of injury or damage to property. The 
% software developed by NIST employees is not subject to copyright 
% protection within the United States.
%
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Refactored code
%              Steve Blandino <steve.blandino@nist.gov>
%              Yunqi Feng <Yunqi.Feng@UGent.be>

%% Input Parameters Management and preallocation
nodePosition = nodeCfgInput.nodePosition;
nPAA_centroids = cellfun(@(x) x.nPAA_centroids, nodeCfgInput.paaInfo );
Mpc = cell(paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.numberOfNodes,...
    max(nPAA_centroids),...
    paraCfgInput.totalNumberOfReflections+4,...
    paraCfgInput.numberOfTimeDivisions+1 );
keepBothQDOutput = strcmp(paraCfgInput.outputFormat, 'both'); 
isJsonOutput = strcmp(paraCfgInput.outputFormat, 'json');
displayProgress = 1;
ts = paraCfgInput.totalTimeDuration/paraCfgInput.numberOfTimeDivisions;
outputPaa = cell(paraCfgInput.numberOfNodes,paraCfgInput.numberOfNodes);
% List of paths
% outputPath = strcat('output/',paraCfgInput.environmentFileName(1:end-4),'_',strrep(datestr(datetime), ':', '_'));
% mkdir(outputPath);
outputPath = strcat('output/',paraCfgInput.environmentFileName(1:end-4),'_');
if exist(outputPath, 'dir')
    [~, ~, ~] =rmdir(outputPath,'s'); 
end
mkdir(outputPath);

ns3Path = fullfile(outputPath, 'Ns3');
qdFilesPath = fullfile(ns3Path, 'QdFiles');
if ~isfolder(qdFilesPath)
    mkdir(qdFilesPath)
end

if paraCfgInput.switchSaveVisualizerFiles == 1
    visualizerPath = fullfile(outputPath, 'Visualizer');
    if ~isfolder(visualizerPath)
        mkdir(visualizerPath)
    end
end

% Init output files
if ~isJsonOutput || keepBothQDOutput
    fids = getQdFilesIds(qdFilesPath, paraCfgInput.numberOfNodes,...
        paraCfgInput.useOptimizedOutputToFile);
end

%% Init
% switchPolarization = 0;
% switchCp = 0;
% polarizationTx = [1, 0];

MaterialLibrary = importMaterialLibrary(paraCfgInput.materialLibraryPath);

% Extracting CAD file and storing in an XMl file, CADFile.xml
[CADop, switchMaterial,scatterer_confg] = getCadOutput(paraCfgInput.environmentFileName,...
    outputPath, MaterialLibrary, paraCfgInput.referencePoint,...
    paraCfgInput.selectPlanesByDist, paraCfgInput.indoorSwitch,...
    paraCfgInput.confg_dim,paraCfgInput.scatterer_num,...
    paraCfgInput.scatterer_dim,paraCfgInput.scatterer_dis, ...
    ref_nodeLocTx,ref_nodeLocRx, ...
    scatterer_confg,paraCfgInput.scatterer_height_max,paraCfgInput.mat_confg,...
    paraCfgInput.antenna,paraCfgInput.Ptx,paraCfgInput.polarization,paraCfgInput.scattering);
CADop_diffr = generate_edges(CADop);
wedges = get_wedge_information(CADop_diffr, CADop);
if paraCfgInput.switchSaveVisualizerFiles == 1
    % Save output file with room coordinates for visualization
    if ~isempty(CADop)
        RoomCoordinates = CADop(:, 1:9);
        csvwrite(fullfile(visualizerPath, 'RoomCoordinates.csv'),...
            RoomCoordinates);
    else
        %RoomCoordinates = CADop(:, 1:9);
        csvwrite(fullfile(visualizerPath, 'RoomCoordinates.csv'),...
            CADop);
    end
end
XPD = 30;
XPD_r = 10^(-XPD/10);
switch(paraCfgInput.polarization)
    case 'V-V'
        Jones = [XPD_r;1];
    case 'H-H'
        Jones = [1;XPD_r];
    case 'dual'
        Jones =[1 XPD_r;XPD_r 1];
end

tmp_plane = [];
plane_mat = [];

%% Loop over time instances
for iterateTimeDivision = 1:paraCfgInput.numberOfTimeDivisions
    if mod(iterateTimeDivision,100)==0 && displayProgress 
        disp([fprintf('%2.2f', iterateTimeDivision/paraCfgInput.numberOfTimeDivisions*100),'%'])
    end
             
    %% Point rotation
    % PAAs not centered [0,0,0] have a
    % different position in the global frame if the node rotates. Compute
    % the new PAAs position as well as the equivalent angle resulting from
    % successive transformations (initial PAA orientation + rotation of the
    % node over time)
    for nodeId = 1:paraCfgInput.numberOfNodes
            centerRotation = nodePosition(iterateTimeDivision,:, nodeId);
            nodeRotationEucAngles = nodeCfgInput.nodeRotation(iterateTimeDivision,:, nodeId);
            paaInitialPosition = reshape(squeeze(...
                nodeCfgInput.paaInfo{nodeId}.centroidTimePosition(iterateTimeDivision,:,:)), [], 3);
            [paaRotatedPosition, nodeEquivalentRotationAngle] = coordinateRotation(paaInitialPosition, ...
                centerRotation,...
                nodeRotationEucAngles ...
                );
            nodeCfgInput.nodeEquivalentRotationAngle(iterateTimeDivision,:, nodeId) = nodeEquivalentRotationAngle;
            nodeCfgInput.paaInfo{nodeId}.centroid_position_rot(iterateTimeDivision,:,:) =paaRotatedPosition;
    end   
    
    %% Iterates through all the PAA centroids
    for iterateTx = 1:paraCfgInput.numberOfNodes
        
        for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
            
            for iteratePaaTx = 1:nPAA_centroids(iterateTx)
                
                for iteratePaaRx = 1:nPAA_centroids(iterateRx)
                    output = [];
                        
                    % Update centroids position
                    Tx = squeeze(nodeCfgInput.paaInfo{iterateTx}.centroid_position_rot(iterateTimeDivision,iteratePaaTx,:)).';
                    Rx = squeeze(nodeCfgInput.paaInfo{iterateRx}.centroid_position_rot(iterateTimeDivision,iteratePaaRx,:)).';
                    
                    % Update rotation Tx struct
                    QTx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateTx);
                    QTx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateTx);
                    
                    % Update rotation Rx struct
                    QRx.center(1,:) = nodePosition(iterateTimeDivision,:,iterateRx);
                    QRx.angle(1,:) = nodeCfgInput.nodeRotation(iterateTimeDivision,:, iterateRx);
                    
                    % Update node velocity
                    previousTxPosition =  squeeze(nodeCfgInput.paaInfo{iterateTx}.centroid_position_rot(max(iterateTimeDivision-1,1),iteratePaaTx,:)).';
                    previousRxPosition =  squeeze(nodeCfgInput.paaInfo{iterateRx}.centroid_position_rot(max(iterateTimeDivision-1,1),iteratePaaRx,:)).';
                    
                    vtx = (Tx-previousTxPosition)./ts;
                    vrx = (Rx-previousRxPosition)./ts;
  
                    % LOS Path generation
                    [isLos, delayLos, output] = LOSOutputGenerator(CADop,CADop_diffr, Rx, Tx,...
                        output, vtx, vrx, paraCfgInput.carrierFrequency, ...
                        paraCfgInput.indoorSwitch,paraCfgInput.enablePhaseOutput, ...
                        paraCfgInput.antenna,paraCfgInput.Ptx, ...
                        paraCfgInput.polarization,Jones,'qTx', QTx, 'qRx', QRx);
                    if isLos
                        LOS_output = output;
                        pseudo_LOS_output = [];
                    else
                        LOS_output=[];
                        pseudo_LOS_output = LOS_output;
                    end
                    if ~isempty(CADop)
                    % Store MPC
                        if paraCfgInput.switchSaveVisualizerFiles && isLos
                            multipath1 = [Tx, Rx];
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, 1, iterateTimeDivision+1} =multipath1;
                        end
                    
                    % Higher order reflections (Non LOS)
                        for iterateOrderOfReflection = 1:paraCfgInput.totalNumberOfReflections
                            numberOfReflections = iterateOrderOfReflection;
                        
                            [ArrayOfPoints, ArrayOfPlanes, numberOfPlanes,...
                                ~, ~, arrayOfMaterials, ~] = treetraversal(CADop,...
                                numberOfReflections, numberOfReflections,...
                                0, 1, 1, 1, Rx, Tx, [], [],...
                                switchMaterial, [], 1);
                        
                            numberOfPlanes = numberOfPlanes - 1;
                        
                            [outputTemporary, multipathTemporary,tmp_plane, plane_mat] = ...
                                multipath(delayLos,...
                                ArrayOfPlanes, ArrayOfPoints, Rx, Tx, ...
                                CADop,CADop_diffr, numberOfPlanes, ...
                                MaterialLibrary, arrayOfMaterials, ...
                                switchMaterial, vtx, vrx, ...
                                paraCfgInput.switchDiffuseComponent,...
                                paraCfgInput.switchQDModel,...
                                paraCfgInput.inputScenarioName(10:end),...
                                paraCfgInput.carrierFrequency,...
                                paraCfgInput.indoorSwitch,...
                                paraCfgInput.enablePhaseOutput,...
                                paraCfgInput.diffusePathGainThreshold,...
                                paraCfgInput.mat_confg,...
                                paraCfgInput.antenna,paraCfgInput.Ptx,...
                                paraCfgInput.polarization,Jones, ...
                                paraCfgInput.scattering,LOS_output,pseudo_LOS_output, ...
                                tmp_plane, plane_mat,...
                                'qTx', QTx, 'qRx', QRx, ...
                                'reflectionLoss', paraCfgInput.reflectionLoss);
                        
                            nMpc = size(multipathTemporary,1);
                            %Store MPC
                            if paraCfgInput.switchSaveVisualizerFiles &&...
                                    nMpc > 0                            
                                    multipath1 = multipathTemporary(:,...
                                        2:end); %Discard reflection order column
                                    Mpc{iterateTx,iteratePaaTx,...
                                        iterateRx,iteratePaaRx, ...
                                        iterateOrderOfReflection+1, iterateTimeDivision+1} =multipath1;                            
                            end
                        
                            %Store QD output                        
                            if size(output) > 0
                                output = [output;outputTemporary]; %#ok<AGROW>
                            
                            elseif size(outputTemporary) > 0
                                output = outputTemporary;
                            
                            end
                            % If no paths found, check for NLOS diffracted paths
                            % if isempty(output)
                            %     [outputTemporary, multipathTemporary] = diffractionPathGenerator(Tx, Rx, CADop_diffr, CADop, vtx, vrx, paraCfgInput.carrierFrequency);
                            %     if ~isempty(outputTemporary)
                            %         output = [output; outputTemporary];
                            %         Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, 1, iterateTimeDivision+1} = multipathTemporary;
                            %     end
                            % end

                        end
                        [outputDiffraction, multipathDiffraction, tmp_plane, plane_mat] = diffractionPathGenerator(Tx, Rx, ...
                            wedges, CADop, vtx, vrx, paraCfgInput.carrierFrequency, ...
                            LOS_output,pseudo_LOS_output,paraCfgInput.polarization,Jones, ...
                            paraCfgInput.Ptx, paraCfgInput.antenna, paraCfgInput.enablePhaseOutput, ...
                            MaterialLibrary, paraCfgInput.mat_confg, paraCfgInput.indoorSwitch, tmp_plane, plane_mat);
                        if paraCfgInput.switchSaveVisualizerFiles && ~isempty(outputDiffraction)
                        % Add the found diffracted paths to the list of multipath components
                            output = [output; outputDiffraction];
                            %outputDiffraction1 = outputDiffraction(:,2:end); %Discard reflection order column
                            % Store the diffracted paths for visualization
                            % NOTE: You'll need to decide how to store this in your Mpc array. 
                            % A simple way is to assign it to a new, higher-order reflection index.
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, paraCfgInput.totalNumberOfReflections+2, iterateTimeDivision+1} = multipathDiffraction;
                        end

                        [outputdifref, multipathdifref, tmp_plane, plane_mat] = ...
                            reflectionDiffractionPathGenerator(Tx, Rx, wedges, CADop, vtx, vrx, ...
                            paraCfgInput.carrierFrequency, LOS_output,pseudo_LOS_output,paraCfgInput.polarization,Jones, ...
                            paraCfgInput.Ptx, paraCfgInput.antenna, paraCfgInput.enablePhaseOutput, ...
                            MaterialLibrary, paraCfgInput.mat_confg, paraCfgInput.indoorSwitch, tmp_plane, plane_mat, ...
                            paraCfgInput.scattering);
                        if paraCfgInput.switchSaveVisualizerFiles && ~isempty(outputdifref)
                        % Add the found diffracted paths to the list of multipath components
                            output = [output; outputdifref];
                            %outputDiffraction1 = outputDiffraction(:,2:end); %Discard reflection order column
                            % Store the diffracted paths for visualization
                            % NOTE: You'll need to decide how to store this in your Mpc array. 
                            % A simple way is to assign it to a new, higher-order reflection index.
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, paraCfgInput.totalNumberOfReflections+3, iterateTimeDivision+1} = multipathdifref;
                        end


                        [outputrefdif, multipathrefdif, tmp_plane, plane_mat] = ...
                            diffractionReflectionPathGenerator(Tx, Rx, wedges, CADop, vtx, vrx, ...
                            paraCfgInput.carrierFrequency, LOS_output,pseudo_LOS_output,paraCfgInput.polarization,Jones, ...
                            paraCfgInput.Ptx, paraCfgInput.antenna, paraCfgInput.enablePhaseOutput, ...
                            MaterialLibrary, paraCfgInput.mat_confg, paraCfgInput.indoorSwitch, tmp_plane, plane_mat, ...
                            paraCfgInput.scattering);
                        if paraCfgInput.switchSaveVisualizerFiles && ~isempty(outputdifref)
                        % Add the found diffracted paths to the list of multipath components
                            output = [output; outputrefdif];
                            %outputDiffraction1 = outputDiffraction(:,2:end); %Discard reflection order column
                            % Store the diffracted paths for visualization
                            % NOTE: You'll need to decide how to store this in your Mpc array. 
                            % A simple way is to assign it to a new, higher-order reflection index.
                            Mpc{iterateTx,iteratePaaTx,iterateRx,iteratePaaRx, paraCfgInput.totalNumberOfReflections+4, iterateTimeDivision+1} = multipathrefdif;
                        end

                        %Tx/Rx antennas will point on the best NLOS direction
                        if strcmp(paraCfgInput.antenna,'dir') && isempty(LOS_output) && ~isempty(output)
                            output_tmp = [];
                            [~,row_ind] = max(output(:,9));
                            AOD_best = output(row_ind,10);
                            EOD_best = output(row_ind,11);
                            AOA_best = output(row_ind,12);
                            EOA_best = output(row_ind,13);
                            [r, ~] = size(output);
                            for ii = 1:r
                                if ii~=row_ind
                                    rel_aaod=abs(output(ii,10)-AOD_best);
                                    rel_eaod=abs(output(ii,11)-EOD_best);
                                    rel_aaoa=abs(output(ii,12)-AOA_best);
                                    rel_eaoa=abs(output(ii,13)-EOA_best);
                                    Gt=horn_gain(rel_eaod,rel_aaod);
                                    Gr=horn_gain(rel_eaoa,rel_aaoa);
                                    output(ii,9)=output(ii,9)-23.7-23.7+Gt+Gr;
                                end
                            end
                            for ii = 1:r
                                if output(ii,9)>-120
                                    output_tmp1 = output(ii,:);
                                    output_tmp = [output_tmp;output_tmp1];
                                end
                            end
                            output = output_tmp;
                        end

                    
                    % Create outputPAA array of struct. Each entry of the
                    % array is a struct relative to a NodeTx-NodeRx 
                    % combination. Each struct has the entries 
                    % - paaTxXXpaaRxYY: channel between paaTx XX and paaRx
                    % YY.  
                        outputPaa{iterateTx, iterateRx}.(sprintf('paaTx%dpaaRx%d', iteratePaaTx-1, iteratePaaRx-1))= output;
                        outputPaa{iterateRx, iterateTx}.(sprintf('paaTx%dpaaRx%d', iteratePaaRx-1, iteratePaaTx-1))= reverseOutputTxRx(output);
                    else
                        outputPaa{iterateTx, iterateRx}.(sprintf('paaTx%dpaaRx%d', iteratePaaTx-1, iteratePaaRx-1))= output;
                        outputPaa{iterateRx, iterateTx}.(sprintf('paaTx%dpaaRx%d', iteratePaaRx-1, iteratePaaTx-1))= reverseOutputTxRx(output);
                    end
                end
            end
        end
    end
    
    %% Generate channel for each PAA given the channel of the centroids
    outputPaaTime(:,:,iterateTimeDivision) = generateChannelPaa(outputPaa, nodeCfgInput.paaInfo);  %#ok<AGROW>
    %outputPaa={[],outputPaa{1,2}.paaTx0paaRx0;outputPaa{2,1}.paaTx0paaRx0,[]};
    
    %% Write QD output in CSV files
    if ~isJsonOutput || keepBothQDOutput
        for iterateTx = 1:paraCfgInput.numberOfNodes
            for iterateRx = iterateTx+1:paraCfgInput.numberOfNodes
                writeQdFileOutput(outputPaaTime{iterateTx, iterateRx,iterateTimeDivision},...
                    paraCfgInput.useOptimizedOutputToFile, fids, iterateTx, iterateRx,...
                    qdFilesPath, paraCfgInput.qdFilesFloatPrecision);
                
                writeQdFileOutput(outputPaaTime{iterateRx,iterateTx,iterateTimeDivision},...
                    paraCfgInput.useOptimizedOutputToFile, fids, iterateRx,iterateTx,...
                    qdFilesPath, paraCfgInput.qdFilesFloatPrecision);
            end
        end
    end
    
    %clear outputPAA
end

%% Write output in JSON files
% QD output
if isJsonOutput  || keepBothQDOutput
    writeQdJsonOutput(outputPaaTime,cellfun(@(x) x.nPaa,  nodeCfgInput.paaInfo),...
        qdFilesPath);
end

if paraCfgInput.switchSaveVisualizerFiles
    Mpc(:,:,:,:,:,1) = [];
    writeVisualizerJsonOutput(visualizerPath, paraCfgInput, nodeCfgInput, nPAA_centroids, nodePosition, Mpc)
end

if ~isJsonOutput || keepBothQDOutput
    closeQdFilesIds(fids, paraCfgInput.useOptimizedOutputToFile);
end

%% Write useful output information. 
writeReportOutput = 0 ; %Set to 0 to allow succeful test.
if writeReportOutput
    f = fopen(strcat(outputPath, filesep,'report.dat'), 'w'); %#ok<UNRCH> 
    fprintf(f, 'Device Rotation:\t%d\n', paraCfgInput.isDeviceRotationOn);
    fprintf(f, 'Initial Orientation:\t%d\n', paraCfgInput.isInitialOrientationOn);
    fprintf(f, 'PAA centered:\t%d\n', paraCfgInput.isPaaCentered);
    fclose(f);
end
end