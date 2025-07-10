function [qdRay, multipath,tmp_plane, plane_mat] =...
    multipath(delayLos, ArrayOfPlanes, ArrayOfPoints, Rx, Tx, CADOutput,...
    CADop_diffr,numberOfRowsArraysOfPlanes, MaterialLibrary, arrayOfMaterials,...
    switchMaterial, velocityTx, velocityRx, ...
    diffuseGeneratorSwitch, qdModelSwitch, scenarioName, frequency,IndoorSwitch, ...
    enablePhase,diffusePathGainThreshold,mat_confg,...
    antenna,Ptx,polarization,Jones,...
    scattering,LOS_output,pseudo_LOS_output,tmp_plane, plane_mat,varargin)
% INPUTS -
% ArrayOfPlanes - Similar to Array of points. Each triangle occupies 4
%   columns (plane equation). The first column has the order of reflection
%   (o/p of treetraversal)
% ArrayOfPoints - combinations of multiple triangles, every row is a unique
%   combination. every triangle occupies 9 columns (3 vertices). (o/p of
%   treetraversal)
% Rx - Rx position
% Tx - Tx position
% CADOutput - CAD output
% numberOfRowsArraysOfPlanes - number of arrays of planes
% MaterialLibrary - Similar to Array of points. Each triangle occupies 1
%   triangle. The data is the row number of material from Material library
%   arrayOfMaterials - Similar to Array of points. Each triangle occupies 1
%   triangle. The data is the row number of material from Material library
% arrayOfMaterials - array of material corresponding to each planes
% switchMaterial - whether triangle materials properties are present
% velocityTx, velocityRx are velocities of tx and rx respectively
% diffuseGeneratorSwitch - Switch to turn ON or OFF the Qausi dterministic 
% module 1 = ON, 0 = OFF
% diffuseGeneratorSwitch - Switch to select the approach to generate   
% diffuse components 
% 2 = approach based on IMEC DBARC model
% 1 = approach based on NIST measured QD parameters, 
% 0 = approach based on channel document QD parameters
% qdModelSwitch - Switch to select the QD model 
% scenarioName - define scenario name 
% frequency - the carrier frequency at which the system operates
% diffusePathGainThreshold - This value is used to filter out diffuse
% components for qdModelSwitch = nistMeasurements 
%
% OUTPUTS -
% qdRay - consists of specular and diffuse multipath parameters
% multipath - consists of specular multipath parameters
%
% The phase information in case of presence of polarization information and 
% is encoded in the Jones vector. In case of absence of polarization, order 
% of reflection is multiplied by pi to give phase shift


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

% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Used MATLAB functions 
%   instead of custom ones, vectorized code, improved access to 
%   MaterialLibrary
% Modified by: Neeraj Varshney <neeraj.varshney@nist.gov>, to calculate RL
%   based on new material libraries given in 802.11ay channel document
% Modified by: Yunqi Feng <Yunqi.Feng@UGent.be> to include the D-band
% measurements 
% September 2023

%% Varargin processing 
p = inputParser;
addParameter(p,'indStoc',1)
addParameter(p,'qTx',struct('center', Tx, 'angle', [0 0 0]))
addParameter(p,'qRx',struct('center', Rx, 'angle', [0 0 0]))
addParameter(p,'reflectionLoss',10);
parse(p, varargin{:});
qTx = p.Results.qTx;
qRx = p.Results.qRx;
rl  = p.Results.reflectionLoss;

%% Init
indexMultipath = 1;
indexOutput = 1;
nVarOut = 22;
sizeArrayOfPlanes = size(ArrayOfPlanes);
multipath = [];
c = getLightSpeed;
wavelength = c / frequency;
outputQd = struct('dRay', cell(sizeArrayOfPlanes(1),1), ...
    'rPreCursor', cell(sizeArrayOfPlanes(1),1), ...
    'rPostCursor', cell(sizeArrayOfPlanes(1),1));
%%
if ~isempty(LOS_output)% with the LOS path
    LOS_TOA=LOS_output(8);
    LOS_P=LOS_output(9);
    LOS_AAOD=LOS_output(10);
    LOS_EAOD=LOS_output(11);
    LOS_AAOA=LOS_output(12);
    LOS_EAOA=LOS_output(13);
    LOS_PHA=LOS_output(18);
    LOS_DOP=LOS_output(20);
end
if ~isempty(pseudo_LOS_output)% with the LOS path
    pseudo_LOS_TOA=pseudo_LOS_output(8);
    pseudo_LOS_P=pseudo_LOS_output(9);
    pseudo_LOS_AAOD=pseudo_LOS_output(10);
    pseudo_LOS_EAOD=pseudo_LOS_output(11);
    pseudo_LOS_AAOA=pseudo_LOS_output(12);
    pseudo_LOS_EAOA=pseudo_LOS_output(13);
    pseudo_LOS_PHA=pseudo_LOS_output(18);
    pseudo_LOS_DOP=pseudo_LOS_output(20);
end
if numberOfRowsArraysOfPlanes>0
    orderOfReflection = ArrayOfPlanes(1,1);
    
    % the for loop iterates through all the rows of ArrayOfPlanes,
    % ArrayOfPoints and provides a single row as input to
    % singleMultipathGenerator function
    % QD model is present in  this loop
    multipath = zeros(numberOfRowsArraysOfPlanes,orderOfReflection * 3 + 1);
    % tmp_plane=[];
    % plane_mat=[];
    for iterateNumberOfRowsArraysOfPlanes = numberOfRowsArraysOfPlanes:-1:1
        indexOrderOfReflection = 1;
        multipath(indexMultipath, (indexOrderOfReflection-1)*3 + 1) = orderOfReflection;
        multipath(indexMultipath, (indexOrderOfReflection-1)*3 + 1 + (1:3)) = Rx;
        Reflected = Rx;
        %Incident_angle=[];
        
        % a single row of ArrayOfPlanes,ArrayOfPoints is fed to
        % singleMultipathGenerator function to know whether a path exists. If a
        % path exists then what are vectors that form the path (stored in
        % multipath parameter)                
        [isMpc,~,dod,doa,multipath,distance,dopplerFactor,...
           ~] = singleMultipathGenerator...
            (iterateNumberOfRowsArraysOfPlanes,orderOfReflection,...
            indexOrderOfReflection,ArrayOfPlanes,...
            ArrayOfPoints,Reflected,Rx,Tx,CADOutput,...
            multipath,indexMultipath,velocityTx,velocityRx);
        Incident_angle = angleOfIncidence(multipath(indexMultipath,:));
        % Apply node rotation
        dod_tmp = dod;
        doa_tmp = doa;
        dod = coordinateRotation(dod,[0 0 0], qTx.angle, 'frame');
        doa = coordinateRotation(doa,[0 0 0], qRx.angle, 'frame');
        
        
        % Compute reflection loss
        if isMpc == 1
            if  switchMaterial == 1  
	            [reflectionLoss,diffuseLoss,scatterPaths,tmp_plane,plane_mat] = getDbarcReflectionLoss(MaterialLibrary,... 
		arrayOfMaterials(iterateNumberOfRowsArraysOfPlanes,:),...
                	    multipath(indexMultipath,:),Incident_angle,mat_confg,...
                        IndoorSwitch,polarization,Jones,scattering, ...
                        ArrayOfPlanes(iterateNumberOfRowsArraysOfPlanes,:), ...
                        tmp_plane,plane_mat); 
             else
                % Assumption: r1 loss at each reflection if
                % material is not present in the CAD file
                reflectionLoss = rl * orderOfReflection;
            end
            
        end
                    
        % Corner case: MPC on the edge of triangles would be considered
        % twice. Check if it has been already stored otherwise discard.
        if isMpc == 1
            for i = 1:indexMultipath - 1
                isMpcNonUnique = 1;
                for j = 1:(orderOfReflection * 3) + 6
                    isMpcNonUnique = isMpcNonUnique && ....
                        (multipath(i,j) == multipath(indexMultipath,j));
                end
                isMpc = isMpc && ~isMpcNonUnique;
            end
        end
        
        % the delay, AoA, AoD, path loss of the path are stored 
        % in output parameter
        dRay = zeros(1, nVarOut);
        if  isMpc == 1
            dRay(1) = indexMultipath;
            % dod - direction of departure
            dRay(2:4) = dod;
            % doa - direction of arrival
            dRay(5:7) = doa;
            % Time delay
            dRay(8) = distance/c;
            
            %dRay(9) = 20*log10(wavelength / (4*pi*distance)) ...
                        %- reflectionLoss;            
            % Aod azimuth
            dRay(10)=atan2(dod(2),dod(1))*180/pi;
            % Aod elevation
            dRay(11)=atan2(dod(3),sqrt(dod(1)^2+dod(2)^2))*180/pi;
            % Aoa azimuth
            dRay(12)=atan2(doa(2),doa(1))*180/pi;
            % Aoa elevation
            dRay(13)=atan2(doa(3),sqrt(doa(1)^2+doa(2)^2))*180/pi;

            % Friis transmission loss
             if ~isempty(LOS_output)% still assume LOS path exists
                 rel_aaod=abs(dRay(10)-LOS_AAOD);
                 rel_eaod=abs(dRay(11)-LOS_EAOD);
                 rel_aaoa=abs(dRay(12)-LOS_AAOA);
                 rel_eaoa=abs(dRay(13)-LOS_EAOA);
                 if strcmp(antenna, 'dir')
                     Gt=horn_gain(rel_eaod,rel_aaod);
                     Gr=horn_gain(rel_eaoa,rel_aaoa);
                 elseif strcmp(antenna, 'omni')
                     Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                     Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
                 elseif strcmp(antenna, 'scan')
                     Gt=23.7;
                     Gr=23.7;
                 end
             elseif isempty(LOS_output) && (strcmp(antenna, 'dir') || strcmp(antenna, 'scan')) 
                 Gt=23.7;
                 Gr=23.7;
 % %If LOS path does not exist in directional channel 
 % sounding schemes, the Tx/Rx antennas will point on the 
 % best NLOS direction, but it's done in post processing
                 % LOS_AAOD=dRay(10);
                 % LOS_EAOD=dRay(11);
                 % LOS_AAOA=dRay(12);
                 % LOS_EAOA=dRay(13);
                 % LOS_output=[LOS_AAOD,LOS_EAOD,LOS_AAOA,LOS_EAOA];
             elseif isempty(LOS_output) && strcmp(antenna, 'omni')
                 rel_aaod=abs(dRay(10)-pseudo_LOS_AAOD);
                 rel_eaod=abs(dRay(11)-pseudo_LOS_EAOD);
                 rel_aaoa=abs(dRay(12)-pseudo_LOS_AAOA);
                 rel_eaoa=abs(dRay(13)-pseudo_LOS_EAOA);
                 Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                 Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
             end

            if strcmp(polarization, 'V-V') || strcmp(polarization, 'H-H')
                dRay(9)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*distance))-reflectionLoss;
                dRay(16)=dRay(9)-Ptx-Gt-Gr;
                Jones_dB = 10*log10(Jones);
                dRay(17) = NaN;%cross-polarized channel gain,reflection
                dRay(19) = min(Jones_dB + dRay(9));
                dRay(22) = NaN;
            %elseif strcmp(polarization, 'cross') || strcmp(polarization, 'dual')
            elseif strcmp(polarization, 'dual')
                %XPD = 34;
                Jones_dB = 10*log10(Jones);
                dRay(9)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*distance))-reflectionLoss(1);%co-polarized channel gain
                dRay(16)=dRay(9)-Ptx-Gt-Gr;
                dRay(17)=Ptx+Gt+Gr+20*log10(wavelength/(4*pi*distance))-reflectionLoss(2);%another co-polarization channel gain
                dRay(19)= min(Jones_dB(:,1) + dRay(9));% cross-polarization channel gain at Tx
                dRay(22)= min(Jones_dB(:,1) + dRay(17)); % cross-polarization channel gain at Rx
            end
            % if dRay(17)<=-120
            %     dRay(17)=NaN;
            % end
            if dRay(19)<=-120
                dRay(19)=NaN;
            end
            if dRay(22)<=-120
                dRay(22)=NaN;
            end
            if enablePhase
                dRay(18) = mod(orderOfReflection*pi+distance/wavelength*2*pi,2*pi);
            else
                dRay(18) = 0;
            end
            dRay(14)=Gt;
            dRay(15)=Gr;
            dRay(20) = dopplerFactor * frequency;
            dRay(21) = 0;
            if dRay(9)<=-120 && (dRay(17)<=-120 || isnan(dRay(17)))
                dRay=[];
                outputQd(indexOutput).dRay=dRay;
                multipath(iterateNumberOfRowsArraysOfPlanes,:)=[];
            else
                if dRay(17)<=-120
                    dRay(17)=NaN;
                end
                if dRay(9)<=-120
                    dRay(9)=NaN;
                end
                outputQd(indexOutput).dRay = dRay;
            end
            % If dRay detectable, generate components
            if  switchMaterial == 1 && scattering == 1 && ~isempty(dRay)
                % [~, rPreCursor, rPostCursor] =...
                %     qdGenerator(delayLos, outputQd(indexOutput).dRay,...
                %     arrayOfMaterials(iterateNumberOfRowsArraysOfPlanes,:),...
                %     MaterialLibrary,qdModelSwitch, scenarioName,...
                %     diffusePathGainThreshold); 
                % 
                % outputQd(indexOutput).rPreCursor  = rPreCursor;
                % outputQd(indexOutput).rPostCursor = rPostCursor;
                % indexOutput = indexOutput + 1;
                % indexMultipath = indexMultipath + 1;
            elseif switchMaterial == 1 && diffuseGeneratorSwitch == 0 && ~isempty(dRay)
                indexOutput = indexOutput + 1;
                indexMultipath = indexMultipath + 1;
            elseif isempty(dRay)
                outputQd(indexOutput).rPreCursor  =[];
                outputQd(indexOutput).rPostCursor =[];
            end
        end
    end
    
    qdRay =     [ ...
    reshape([outputQd.dRay],nVarOut, []).';...
    reshape([outputQd.rPreCursor].', nVarOut, []).';...
    reshape([outputQd.rPostCursor].', nVarOut, []).'];
    qdRay(isnan(qdRay(:,1)),:) =[];
    
    if indexMultipath>=1
        multipath(indexMultipath:end,:) = [];
    end
    
end
end
function incidentAngle = angleOfIncidence(multipath) 
% ANGLEOFINCIDENCE returns angle of incident for first and second order
% reflections
% 
% Inputs:
% multipath - consists of specular multipath parameters. This vector is
% used to calculate angle(s) of incident.
% 
% Output: 
% incidentAngle - incident angles till second order reflections

differenceVectorRxFirstPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
% differenceVectorRxFirstPoI is the difference vector between Rx 
% and first point of incidence (PoI).
differenceVectorTxFirstPoI = (multipath(1,8:10) - multipath(1,5:7))...
                            /( norm(multipath(1,8:10) - multipath(1,5:7)));
% differenceVectorTxFirstPoI is the difference vector between Tx 
% and first point of incidence (PoI). 
dpAoI = dot(differenceVectorRxFirstPoI, differenceVectorTxFirstPoI);
% dpAoI is the dot product between differenceVectorRxFirstPoI and 
% differenceVectorTxFirstPoI will give the cos of angle between the vectors.
incidentAngle(1) = 0.5*acosd(dpAoI);
% Half of the angle between the vectors differenceVectorRxFirstPoI and 
% differenceVectorTxFirstPoI is the angle of incidence. This is because 
% angle of incidence is equal to angle of reflection.
if multipath(1,1) == 2 % This is for second order reflection
    differenceVectorRxSecondPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
    differenceVectorFirstPoISecondPoI =...  
                            (multipath(1,8:10) - multipath(1,5:7))...
                            /( norm(multipath(1,8:10) - multipath(1,5:7)));
    differenceVectorSecondPoIFirstPoI = -differenceVectorFirstPoISecondPoI;
    differenceVectorTxFirstPoI =...
        (multipath(1,11:13) - multipath(1,8:10))/...
        norm((multipath(1,11:13) - multipath(1,8:10)));
    dpAoI = dot(differenceVectorSecondPoIFirstPoI,...
                differenceVectorTxFirstPoI);
    incidentAngle(1) = 0.5*acosd(dpAoI);
    dpAoI = dot(differenceVectorRxSecondPoI,...
                differenceVectorFirstPoISecondPoI);
    incidentAngle(2)  = 0.5*acosd(dpAoI);
end
if multipath(1,1) == 3 % This is for third order reflection
    differenceVectorRxThirdPoI = (multipath(1,2:4) - multipath(1,5:7))...
                                /norm(multipath(1,2:4) - multipath(1,5:7));
    differenceVectorRxSecondPoI = (multipath(1,5:7) - multipath(1,8:10))...
                            /( norm(multipath(1,5:7) - multipath(1,8:10)));
    differenceVectorFirstPoISecondPoI = (multipath(1,8:10) - multipath(1,11:13))...
                            /( norm(multipath(1,8:10) - multipath(1,11:13)));
    differenceVectorTxFirstPoI =...
        -(multipath(1,11:13) - multipath(1,14:16))/...
        norm((multipath(1,11:13) - multipath(1,14:16)));
    dpAoI = dot(differenceVectorTxFirstPoI,...
                differenceVectorFirstPoISecondPoI);
    incidentAngle(1) = 0.5*acosd(dpAoI);
    dpAoI = dot(-differenceVectorFirstPoISecondPoI,...
                differenceVectorRxSecondPoI);
    incidentAngle(2)  = 0.5*acosd(dpAoI);
    dpAoI = dot(-differenceVectorRxSecondPoI,...
                differenceVectorRxThirdPoI);
    incidentAngle(3)  = 0.5*acosd(dpAoI);

end
end
