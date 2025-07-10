function [isLOS, delayLOS, output, varargout] = LOSOutputGenerator(CADoutput,CADop_diffr, Rx, Tx,...
    output, velocityTx, velocityRx, frequency, IndoorSwitch,enablePhase,antenna, ...
    Ptx,polarization,Jones,varargin)
% This part of code compute LOS between two nodes
%
% Inputs:
% CADoutput - CAD output
% Tx and Rx locations if using two nodes
% output - multipath parameters
% velocityTx, velocityRx are velocities of tx and rx respectively
% isPolarization - a boolean to describe whether polarization is
%   selected
% isXPol - a boolean to describe whether cross polarization is selected
%   or not. 1 means there is cross polarization and 0 means there is no 
%   cross polarization
% PolarizationTx - gives polarization information of Tx location
% frequency: the carrier frequency at which the system operates
% IndoorSwitch: 0=outdoor,1=indoor
%
%Outputs:
% isLOS - a boolean which gives information whether LOS exist or not.
%   1 stands for existant while 0 is for non existant case
% output - multipath parameters

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
% Modified by: Mattia Lecci <leccimat@dei.unipd.it>, Used MATLAB functions 
%   instead of custom ones
% 2020 NIST/CTL (steve.blandino@nist.gov) 
% Modified by: Yunqi Feng <Yunqi.Feng@UGent.be> to include the D-band
% measurements 
% September 2023

%% Input processing
p = inputParser;
addParameter(p,'qTx',struct('center', Tx, 'angle', [0 0 0]))
addParameter(p,'qRx',struct('center', Rx, 'angle', [0 0 0]))
parse(p, varargin{:});
qTx = p.Results.qTx;
qRx = p.Results.qRx;
% Direction of departure (DoD) is simple the difference of position vectors
% of Tx and Rx
dodNoRot = Rx - Tx;
dod = coordinateRotation(dodNoRot, [0 0 0], qTx.angle, 'frame');
% distance is the total length of multipath
distance=norm(dod);
% delayLOS is the delay of LOS path
c=getLightSpeed;
delayLOS = distance/c;
% Direction of arrival (DoA) is negative of DoD
%doaNoRot = Tx - Rx;
doaNoRot = -1 * (Tx - Rx);
doa = coordinateRotation(doaNoRot, [0 0 0], qRx.angle,'frame');
% Calculating Doppler factor for LOS
velocityTxAlongDirectionOfDeparture=dot(velocityTx,-1.*dod);
velocityRxAlongDirectionOfDeparture=dot(velocityRx,-1.*dod);
dopplerFactor=(velocityRxAlongDirectionOfDeparture...
    -velocityTxAlongDirectionOfDeparture)/c;
% To verify whether DoA vector exists
if IndoorSwitch == 0
    std = 2;
elseif IndoorSwitch == 1
    std = 1;
end
SF = normrnd(0,std);
isLOS = verifyPath(Tx, Rx, doaNoRot, [0,0,0],...
    [0,0,0], CADoutput, 2, false);

    output1 = nan(1,22);
    
    lambda=c/frequency;
    output1(1) = 1;
    % dod - direction of departure
    output1(2:4) = dod;
    % doa - direction of arrival
    output1(5:7) = doa;
    % Time delay
    output1(1,8) = delayLOS;
    % In LOS path, assume antennas are exactly pointing each other
    if strcmp(antenna, 'omni')
        Gt=3+rand;
        Gr=3+rand;
    else
        Gt=23.7;
        Gr=23.7;
    end
    if isLOS==1 % if DoA exists
        % Path gain
        if strcmp(polarization, 'V-V') || strcmp(polarization, 'H-H')
            Jones_dB = 10*log10(Jones);
            output1(1,9) = Ptx+Gt+Gr+20*log10(lambda/(4*pi*distance)) + SF;%co-polarized channel gain          
            output1(1,16) = output1(1,9)-(Ptx+Gt+Gr);%path loss
            output1(1,17) = NaN;%another co-polarized channel gain
            output1(1,19) = min(Jones_dB + output1(1,9));%cross-polarized channel gain
            output1(1,22) = NaN;
        elseif strcmp(polarization, 'dual')
            Jones_dB = 10*log10(Jones);
            output1(1,9) = Ptx+Gt+Gr+20*log10(lambda/(4*pi*distance));
            output1(1,16) = output1(1,9)-(Ptx+Gt+Gr);
            output1(1,17) = output1(1,9);% dual-polarized channel gain
            output1(1,19) = min(Jones_dB(:,1) + output1(1,9));
            output1(1,22) = min(Jones_dB(:,1) + output1(1,9));
        end
        if output1(9)<=-120
            isLOS=0;
        end
        if output1(1,17)<=-120
            output1(1,17) =NaN;
        end
        if output1(1,19)<-120
            output1(1,19)=NaN;
        end
        if output1(1,22)<-120
            output1(1,22)=NaN;
        end
        [diffraction_loss, diffraction_points] = calculate_diffraction_on_path(Tx, Rx, CADop_diffr, lambda);

    % Add diffraction loss to the path gain
        output1(1,9) = output1(1,9) - diffraction_loss;

    % Store diffraction points for visualization (optional)
        if ~isempty(diffraction_points)
            varargout{1} = diffraction_points;
        end
    else
        %warning('No direct path exists, adding OBS path with PL 150 dB');
        %output1(1,9) = -150;
        output1(9)=-120;
    end
    % Aod azimuth
    output1(10)=atan2(dod(2),dod(1))*180/pi;
    % Aod elevation
    output1(11)=atan2(dod(3),sqrt(dod(1)^2+dod(2)^2))*180/pi;
    % Aoa azimuth
    output1(12)=atan2(doa(2),doa(1))*180/pi;
    % Aoa elevation
    output1(13)=atan2(doa(3),sqrt(doa(1)^2+doa(2)^2))*180/pi;
    % Polarization Jones vector
    % if isPolarization
    %     output1(14:15) = PolarizationTx(1,:);
    %     % Cross polarization Jones vector
    %     if isXPol
    %         output1(16:17) = PolarizationTx(2,:);
    %     end
    % end
    output1(14)=Gt;
    output1(15)=Gr;
    if enablePhase
        output1(18) = mod(distance/lambda*2*pi,2*pi);
    else
        output1(18) = 0;
    end
    % Doppler Factor
    output1(20) = dopplerFactor * frequency;
    output1(21) = 0;

    if size(output)>0
        if output1(9)<=-120
            output=output;
        else
            output = [output; output1];
        end
    else
        if output1(9)<=-120
            output=[];
        else
            output = output1;
        end
    end
    
end
