function [scatterPaths] = computeDiffuseScattering(diffuseLoss,normal,doa,dod,materialParams, nRaysDiffuse)
%COMPUTEDIFFUSESCATTERING  Example of Lambertian diffuse scattering (802.11ay style)
%
%   posHit         = [x; y; z] coordinates where ray hits the surface
%   n              = surface normal at posHit (3x1 vector)
%   posTx          = transmitter position (3x1)
%   posRx          = receiver position (3x1) [optional; not always needed here]
%   freqGHz        = frequency in GHz (e.g., 60)
%   materialParams = struct containing .GammaDiff and .specularRefCoeff, etc.
%   nRaysDiffuse   = how many random rays to sample for diffuse
%
%   scatterPaths   = struct array with fields:
%       .direction    (3x1 unit vector for scattered ray)
%       .gain         (power gain or complex amplitude)
%       .phase        (random phase)
%       .delay        (propagation delay from posHit to Rx)
%
%   This function returns 'nRaysDiffuse' random directions following a Lambertian
%   distribution around the surface normal. Each path includes a basic power/delay
%   term. The code is an illustrative *simplification* of 802.11ay-like diffuse
%   modeling. Actual standardization includes further detail and calibration.

    % Unpack some material parameters
    GammaDiff = materialParams;  % Diffuse scattering coefficient
    % specRefCoeff = materialParams.specularRefCoeff;  (if needed for other routines)

    c = 3e8;  % speed of light
    lambda = c / (140 * 1e9);

    % ---------------------------------------------------------
    % 1. Compute incident angle, incident power, etc.
    % ---------------------------------------------------------
    % Incident direction (unit) from Tx to the hit point:
    incidentVec = doa;
    distInc     = norm(incidentVec);
    incDir      = incidentVec / distInc;

    % If you have the power at Tx, Pt, and pathloss up to posHit:
    %   Pr_incident = Pt * somePathLoss
    % For simplicity here, assume we have an incident power of 1.0 at posHit:
    %P_incident = 1.0;

    % Dot with surface normal to find how "oblique" the incidence is
    % [s1,~] = size(normal);
    % if s1 == 1
    %     dotInc = dot(normal, incDir);
    %     if dotInc < 0
    %         dotInc = abs(dotInc); % Normal might be "inside" the surface
    %         normal = -normal;               % (optional) flip if needed
    %     end
    % elseif s1 == 2
    % 
    % dotInc1 = dot(normal1, incDir);
    % if dotInc1 < 0
    %     dotInc1 = abs(dotInc1); % Normal might be "inside" the surface
    %     normal1 = -normal1;               % (optional) flip if needed
    % end
    % dotInc2 = dot(normal2, incDir);
    % if dotInc2 < 0
    %     dotInc2 = abs(dotInc2); % Normal might be "inside" the surface
    %     normal2 = -normal2;               % (optional) flip if needed
    % end
    dotInc = dot(normal, incDir);
    if dotInc < 0
        dotInc = abs(dotInc); % Normal might be "inside" the surface
        normal = -normal;               % (optional) flip if needed
    end

    % ---------------------------------------------------------
    % 2. Generate random outgoing directions for diffuse paths
    %    with Lambertian distribution about the normal n
    % ---------------------------------------------------------
    scatterPaths = struct([]);
    for iRay = 1:nRaysDiffuse

        % Sample random direction in the hemisphere above 'n'
        % Using a cosine-weighted (Lambertian) distribution:
        outDir = sampleLambertianDirection(normal);

        % Check n·outDir
        dotOut = max(dot(normal, outDir), 0);

        % -----------------------------------------------------
        % 3. Compute the fraction of power scattered in outDir
        % -----------------------------------------------------
        % Basic Lambertian:  (GammaDiff * P_incident * (n·incDir)/π ) * (n·outDir)
        % But we distribute over 2π sr, so the factor 1/π is typical for the
        % BRDF, while the hemisphere integration is 2π. The net normalization
        % can vary in different references; here is a common approach:
        diffuseLoss = 10.^(-diffuseLoss./10);
        P_thisRay = diffuseLoss * dotInc * (1/pi) * dotOut;

        % -----------------------------------------------------
        % 4. Compute path gain (linear) or amplitude
        % -----------------------------------------------------
        % Let's treat P_thisRay as *power* for the moment
        % Then amplitude ~ sqrt(P_thisRay)
        amp = sqrt(P_thisRay);

        % Random phase from 0..2π for diffuse
        phs = 2*pi*rand(1);

        % -----------------------------------------------------
        % 5. Compute path to receiver (if needed)
        % -----------------------------------------------------
        % You might send each diffuse ray to the receiver, or
        % just store them for a subsequent geometry check.
        toRxVec = dod;
        distOut = norm(toRxVec);
        % Approximating we keep the same outDir for direction:
        % In practice, you'd see if the outDir actually "points" to the receiver.
        % If you want a strict geometric intersection, you can check the angle
        % between outDir and (posRx - posHit).
        cosAngleRx = dot(outDir, toRxVec/distOut);
        if cosAngleRx < 0
            % This diffuse path does not reach Rx directly
            % (or it hits from the backside). You might skip it or continue
            continue;
        end

        % Additional free-space pathloss from posHit to posRx
        % FSPL(d) = (4πd / λ)^2  =>  power ratio = 1 / FSPL
        fspl = (4*pi*distOut / lambda)^2;
        P_rx = P_thisRay / fspl;
        amp_rx = sqrt(P_rx);

        % Let's combine amplitude with random phase:
        ampTotal = amp_rx;  
        phaseTotal = phs;   % could also add more geometric phase

        % Delay = (distInc + distOut) / c
        delayTotal = (distInc + distOut) / c;

        % Fill out the struct
        scatterPaths(end+1).direction = outDir;
        scatterPaths(end).gain        = 10*log10(ampTotal);
        scatterPaths(end).phase       = phaseTotal;
        scatterPaths(end).delay       = delayTotal;
    end
end

function dirLambert = sampleLambertianDirection(n)
%SAMPLELAMBERTIANDIRECTION  Return a random direction in the hemisphere about n
%
%  Using cosine-weighted hemisphere sampling in local coords:
%      phi   ~ U(0, 2π)
%      r     ~ sqrt(U(0,1))
%      z     ~ sqrt(1 - r^2)
%  Then transform from local (z=normal) to global coordinates.

    % Randoms
    phi   = 2*pi*rand;
    r     = sqrt(rand);
    z     = sqrt(1 - r^2);

    xLocal = r * cos(phi);
    yLocal = r * sin(phi);
    zLocal = z;

    % Build a local coordinate frame: n as z-axis
    % Any orthonormal basis can do; here’s one approach:
    w = n / norm(n);  % w = new z
    % create orthonormal u, v with w
    if abs(w(3)) < 0.999  % if n not close to z-axis
        u = cross([0; 0; 1], w);
        u = u / norm(u);
    else
        u = cross([1; 0; 0], w);
        u = u / norm(u);
    end
    v = cross(w, u);

    % Transform local coords to global
    dirLambert = xLocal*u + yLocal*v + zLocal*w;
    dirLambert = dirLambert / norm(dirLambert);
end
