function demoDiffuseScattering()
% DEMODIFFUSESCATTERING: Example script to show how to incorporate
% diffuse scattering given material props and angle of incidence.

%% 0) Define input parameters
eps_r       = 5.0;        % Relative permittivity
sigma       = 0.01;       % Conductivity (S/m)
freq        = 28e9;       % 28 GHz
thetaIncDeg = 30;         % Incidence angle from normal
polarization = 'TE';      % or 'TM'

% Scattering model parameters
alpha = 0.1;    % fraction of reflected power that goes diffuse
Ndiffuse = 30;  % number of diffuse rays to spawn

%% 1) Compute Fresnel Reflection Coefficients
[Gamma_TE, Gamma_TM] = fresnelCoeffs(eps_r, sigma, freq, thetaIncDeg);

% Pick whichever polarization's reflection you need
switch lower(polarization)
    case 'te'
        R_spec = abs(Gamma_TE).^2;   % power reflectivity
    case 'tm'
        R_spec = abs(Gamma_TM).^2;
    otherwise
        error('Unknown polarization. Choose TE or TM.');
end

fprintf('Specular reflection power coefficient: %.3f\n', R_spec);

%% 2) Partition the power into specular and diffuse
[R_specular_new, R_diffuse] = partitionDiffusePower(R_spec, alpha);
fprintf('  -> New specular portion = %.3f\n', R_specular_new);
fprintf('  -> Diffuse portion      = %.3f\n', R_diffuse);

%% 3) Generate diffuse rays (Lambertian from the surface)
directions = sampleLambertianDirections(Ndiffuse);

% Each diffuse ray has power: R_diffuse / Ndiffuse (discrete approximation)
rayPower = (R_diffuse / Ndiffuse) * ones(Ndiffuse,1);

%% 4) Display results
disp('Diffuse Ray Directions (each row is [x, y, z]):');
disp(directions);
disp('Diffuse Ray Powers:');
disp(rayPower);

% The "specular" ray would be reflected at the angle of reflection = angleInc
%   If the surface normal is +Z and the incidence angle is wrt normal,
%   the specular reflection direction in 3D might be computed from the 
%   incident direction & normal.

% For demonstration, let's just show that we reduced the specular power:
fprintf('\nSpecular ray power is now: %.3f\n', R_specular_new);

end

function directions = sampleLambertianDirections(N)
% SAMPLES LAMBERTIAN directions around the Z-axis (0,0,1).
% Returns an Nx3 matrix of unit direction vectors.

directions = zeros(N,3);

for i = 1:N
    % Random azimuth in [0, 2*pi)
    phi = 2*pi*rand;
    % Random 'r' in [0,1) to determine the polar angle for a cos distribution
    %   We can use the fact that the PDF for Lambertian w.r.t. angle is cos(theta).
    %   A common approach: pick z = sqrt(rand) for half-space distribution.
    z = sqrt(rand); 
    r = sqrt(1 - z^2);
    x = r * cos(phi);
    y = r * sin(phi);
    % z is the "cos(theta)" part
    directions(i,:) = [x, y, z];  % Points in the upper hemisphere
end

% directions are automatically unit vectors in the upper hemisphere
end

function [R_spec_new, R_diffuse] = partitionDiffusePower(R_spec, alpha)
% PARTITIONDIFFUSEPOWER: Splits the reflected power into specular and diffuse
%   R_spec   - specular reflection power coefficient (0..1)
%   alpha    - fraction that goes into diffuse (0..1)

% Ensure alpha is between 0 and 1
alpha = max(0, min(alpha, 1));

% Partition power
R_diffuse  = alpha * R_spec;
R_spec_new = (1 - alpha) * R_spec;

end

function [Gamma_TE, Gamma_TM] = fresnelCoeffs(eps_r, sigma, freq, thetaIncDeg)
% FRESNELCOEFFS: Computes TE/TM reflection coefficients at a planar boundary
%   eps_r       - relative permittivity (real part)
%   sigma       - conductivity (S/m)
%   freq        - frequency (Hz)
%   thetaIncDeg - angle of incidence in degrees (from normal)

% Constants
c0    = 3e8;               % Speed of light in vacuum (m/s)
eps0  = 8.854187817e-12;   % Permittivity of free space
omega = 2*pi*freq;         % Angular frequency
lambda = c0 / freq;

% Complex relative permittivity
eps_complex = eps_r - 1i*(sigma/(omega*eps0));

% Convert angle to radians
theta_inc = deg2rad(thetaIncDeg);

% Snell's law (for the transmitted angle)
% air -> material
n1 = 1;  % air
n2 = sqrt(eps_complex); 
theta_t = asin(n1 .* sin(theta_inc) ./ real(n2));  % approximate real(n2) for transmission

% TE (perpendicular) reflection
% Gamma_TE = (Z2*cos(theta_inc) - Z1*cos(theta_t)) / (Z2*cos(theta_inc) + Z1*cos(theta_t))
% where Z1=eta0, Z2=eta0 / sqrt(eps_complex). We'll do a direct formula.
eta0 = 120*pi;                         % free space impedance ~ 377 ohms
eta2 = eta0 ./ sqrt(eps_complex);     % wave impedance in the material
Gamma_TE = (eta2.*cos(theta_inc) - eta0.*cos(theta_t)) ...
         ./ (eta2.*cos(theta_inc) + eta0.*cos(theta_t));

% TM (parallel) reflection
% Gamma_TM = (eta0*cos(theta_inc) - eta2*cos(theta_t)) / (eta0*cos(theta_inc) + eta2*cos(theta_t))
Gamma_TM = (eta0.*cos(theta_t) - eta2.*cos(theta_inc)) ...
         ./ (eta0.*cos(theta_t) + eta2.*cos(theta_inc));

end



