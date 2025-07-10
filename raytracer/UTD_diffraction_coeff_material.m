function [D_s, D_h] = UTD_diffraction_coeff_material(phi, phi_n, n, L, k, ...
    wedge_face1, wedge_face2, materialLibrary, freq, polarization)
% UTD_DIFFRACTION_COEFF_MATERIAL Calculates UTD coefficients for a wedge with real materials.
%
% This function calculates the Fresnel reflection coefficients for each face of the
% wedge and uses them to modify the standard UTD diffraction coefficients.
%
% Inputs:
%   phi, phi_n, n, L, k - Standard UTD parameters.
%   wedge_face1, wedge_face2 - Structs containing info about each wedge face,
%                              including plane equation and material index.
%   materialLibrary     - The material library.
%   freq                - The frequency of the wave.
%   polarization        - 'V-V' (TE) or 'H-H' (TM).
%
% Outputs:
%   D_s, D_h - Soft and Hard polarization diffraction coefficients.

    % --- 1. Calculate Fresnel Reflection Coefficients for each face ---
    
    % Reflection from face 0 (the face the incident ray hits)
    [R0_s, R0_h] = calculate_fresnel(wedge_face1.mat_idx, phi_n, materialLibrary, freq);

    % Reflection from face n (the other face of the wedge)
    [Rn_s, Rn_h] = calculate_fresnel(wedge_face2.mat_idx, n*pi - phi, materialLibrary, freq);

    % --- 2. Calculate the standard UTD terms ---
    exp_term = exp(1j * (pi / 4));
    sqrt_term = sqrt(8 * pi * k * L);

    beta = [(phi - phi_n), (phi + phi_n), (phi + phi_n) - 2*pi*n, (phi - phi_n) + 2*pi*n];
    % a = 2 * k * L * (cos((2 * pi * n * [0, 1, 2, 3] - beta) / 2)).^2;
    

    F_x = zeros(1, 4); % Pre-allocate
    % F_x = 2j * sqrt(a) .* exp(1j * a) .* integral(@(t) exp(-1j*t.^2), 0, sqrt(a));
    for i = 1:4
        a = 2 * k * L * (sin((beta(i)) / 2)).^2;
        sqrt_a = sqrt(a);
        if sqrt_a > 0
             fresnel_integral_val = integral(@(t) exp(-1j*t.^2), 0, sqrt_a);
             F_x(i) = 2j * sqrt_a * exp(1j * a) * fresnel_integral_val;
        else
            F_x(i) = 0;
        end
    end
    
    % --- 3. Combine with Reflection Coefficients ---

    % For soft polarization (TE - Perpendicular)
    % cot_term_s = cot((pi + beta) / (2 * n));
    % terms_s = [cot_term_s(1) - cot_term_s(2), ...
    %            cot_term_s(3) * Rn_s - cot_term_s(4) * R0_s];
    % D_s = -exp_term / (2 * n * sqrt_term) * sum(terms_s .* F_x);
    % 
    % % For hard polarization (TM - Parallel)
    % cot_term_h = cot((pi - beta) / (2 * n));
    % terms_h = [cot_term_h(1) + cot_term_h(2), ...
    %            cot_term_h(3) * Rn_h + cot_term_h(4) * R0_h];
    % D_h = -exp_term / (2 * n * sqrt_term) * sum(terms_h .* F_x);
    cot_term = cot((pi + beta) / (2 * n));
    sum_term_s = (cot_term(1) * F_x(1)) + (R0_s * cot_term(2) * F_x(2)) + ...
                 (Rn_s * cot_term(3) * F_x(3)) + (R0_s * Rn_s * cot_term(4) * F_x(4));
    D_s = -1 / (2 * n * sqrt(2*pi*k)) * sum_term_s;

    % For hard polarization (TM - Parallel)
    cot_term = cot((pi - beta) / (2 * n)); % Note the sign change in beta for TM
    sum_term_h = (cot_term(1) * F_x(1)) + (R0_h * cot_term(2) * F_x(2)) + ...
                 (Rn_h * cot_term(3) * F_x(3)) + (R0_h * Rn_h * cot_term(4) * F_x(4));
    D_h = -1 / (2 * n * sqrt(2*pi*k)) * sum_term_h;

end

function [Gamma_s, Gamma_h] = calculate_fresnel(mat_idx, incident_angle_rad, materialLibrary, freq)
% CALCULATE_FRESNEL Calculates Fresnel reflection coefficients for a material.
    
    % Get material properties from the library
    permittivity = materialLibrary.permittivity(mat_idx);
    conductivity = materialLibrary.conductivity(mat_idx);
    
    % Constants
    e0 = 8.854e-12;
    mu0 = 4 * pi * 1e-7;
    omega = 2 * pi * freq;
    
    % Complex permittivity
    eps_complex = permittivity * e0 - 1j * conductivity / omega;
    
    % Wave impedances
    eta0 = sqrt(mu0 / e0); % Free space
    eta1 = sqrt(1j * omega * mu0 / (1j * omega * eps_complex)); % Material
    
    % Transmission angle (Snell's Law)
    theta_t = asin(sqrt(e0 * mu0 / real(eps_complex * mu0)) * sin(incident_angle_rad));
    
    % Fresnel equations
    Gamma_h = (eta1 * cos(incident_angle_rad) - eta0 * cos(theta_t)) / ...
              (eta1 * cos(incident_angle_rad) + eta0 * cos(theta_t)); % TM (Parallel)
    Gamma_s = (eta0 * cos(incident_angle_rad) - eta1 * cos(theta_t)) / ...
              (eta0 * cos(incident_angle_rad) + eta1 * cos(theta_t)); % TE (Perpendicular)
end