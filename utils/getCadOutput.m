function [CADop, switchMaterial,scatterer_confg] = getCadOutput(environmentFileName,...
    inputPath, MaterialLibrary, referencePoint, selectPlanesByDist, ...
    indoorSwitch,confg_dim,scatterer_num,scatterer_dim,scatterer_dis, ...
    ref_nodeLocTx,ref_nodeLocRx,scatterer_confg, ...
    scatterer_height_max,mat_confg,antenna,Ptx,polarization,scattering)
%GETCADOUTPUT Function to handle smart CAD file import. It tries to create
% a .mat cache file containing the preprocessed CAD file, as importing a
% raw CAD can be vary time consuming.
% The code checks whether a cache already exist and, if so, whether it is
% older than the CAD file. If the cache does not exist or is outdated,
% import the CAD from scratch a create (or overwrite) the cached output,
% otherwise quickly import the cache.
%
% INPUTS:
% - environmentFileName: file name of CAD environment
% - inputPath: input path of selected scenario, containing environmentFileName
% - MaterialLibrary: see XMLREADER
% - referencePoint: see XMLREADER
% - selectPlanesByDist: see XMLREADER
% - indoorSwitch: see XMLREADER
% - confg_dim: set flexible configuration dimension
%
% OUTPUTS: Same outputs as XMLREADER
%
% SEE ALSO: XMLREADER


% Copyright (c) 2019, University of Padova, Department of Information
% Engineering, SIGNET lab.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


% Importing CAD/AMF/XML files is very slow, especially for large scenarios
% Check first if cached
% Modified by Yunqi Feng <Yunqi.Feng@UGent.be> to return the scatterer
% configuration.
NodeLocTx=ref_nodeLocTx;
NodeLocRx=ref_nodeLocRx;
cacheFilePath = fullfile(inputPath, 'cachedCadOutput.mat');
environmentFilePath = fullfile('environments', environmentFileName);

if exist(cacheFilePath, 'file')
    cacheAttribs = dir(cacheFilePath);
    envirnomentAttribs = dir(environmentFilePath);
    
    % If the cache is older than the env. file, it might have been changed
    % Load it only if cache is recent
    if cacheAttribs.datenum >= envirnomentAttribs.datenum
        load(cacheFilePath,...
            'CADop', 'switchMaterial');
        return
        
    else
        warning(['Cache file exists but might be outdated. ',...
            'Ignoring it and generating new one.']);
        
    end
    
end

% Cache either does not exist or outdated
tmpXmlFilePath = fullfile(inputPath, 'CADFile.xml');

[~,~,extension] = fileparts(environmentFileName);
switch(extension(2:end))
    case {'xml', 'amf'}
        copyfile(environmentFilePath,tmpXmlFilePath);
        
        [CADop, switchMaterial,scatterer_confg] = xmlreader(tmpXmlFilePath,...
            MaterialLibrary, referencePoint, selectPlanesByDist,...
            indoorSwitch,confg_dim,scatterer_num,scatterer_dim, ...
            scatterer_dis,NodeLocTx,NodeLocRx, ...
            scatterer_confg,scatterer_height_max,mat_confg);
        
    case 'obj'
        [CADop, switchMaterial] = importObjFile(environmentFilePath,...
            MaterialLibrary, referencePoint, selectPlanesByDist);
    otherwise
        error('Cannot handle ''%s'' extension properly', extension)
        
end

save(cacheFilePath, 'CADop', 'switchMaterial');

end