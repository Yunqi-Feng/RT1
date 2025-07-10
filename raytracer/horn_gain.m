function G = horn_gain(elevation, azimuth)
% Load the gain matrix from the .mat file
currentDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(currentDir); % Get the parent directory
% Construct the full path to the data file
dataFile = fullfile(parentDir, 'horn_gain.mat'); % Replace with the actual filename
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
