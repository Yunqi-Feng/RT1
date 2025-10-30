function sigma_0 = calculateBeckmannKirchhoff(f, n_complex, sigma, L, theta1_deg, theta2_deg, theta3_deg, polarization, lx, ly)
    % calculateBeckmannKirchhoff Models the Beckmann-Kirchhoff scattering coefficient
    %
    % This function calculates the mean scattered power coefficient <pp*>
    % from a rough surface based on the extended Kirchhoff model presented
    % by Beckmann, as described in the paper "Diffuse Scattering From 
    % Rough Surfaces in THz Communication Channels" (Jansen et al., 2011).
    %
    % Inputs:
    %   f           - Frequency of the incident wave (Hz)
    %   n_complex   - Complex refractive index of the material (e.g., n - 1i*kappa) [cite: 32]
    %   sigma       - Standard deviation of surface height (RMS height) (m) [cite: 69]
    %   L           - Correlation length of the surface (m) [cite: 69]
    %   theta1_deg  - Incident angle (degrees), measured from the normal [cite: 40]
    %   theta2_deg  - Scattering angle (degrees), measured from the normal [cite: 40]
    %   theta3_deg  - Azimuthal scattering angle (degrees) [cite: 40, 50]
    %   polarization - 'TE' or 'TM' (Transverse Electric or Transverse Magnetic) 
    %   lx          - Lateral dimension of the illuminated area (x-direction) (m) [cite: 63]
    %   ly          - Lateral dimension of the illuminated area (y-direction) (m) [cite: 63]
    %
    % Output:
    %   sigma_0     - Mean scattering coefficient <pp*>, a unitless power ratio [cite: 92]
    %
    % Reference Equations (Jansen et al., 2011):
    %   Eq. (3)   : F (Geometric factor) [cite: 65]
    %   Eq. (4)   : <pp*>_inf (Infinitely conductive surface) 
    %   Eq. (5)   : rho_0 (Specular sinc function) [cite: 76]
    %   Eq. (6)-(7): vx, vy (Scattering vector components) [cite: 79, 80]
    %   Eq. (8)   : g (Rayleigh roughness factor) [cite: 81]
    %   Eq. (9)   : r (Fresnel reflection coefficient) [cite: 89, 91]
    %   Eq. (10)  : <pp*> (Finitely conducting surface) [cite: 93]
    
    %% 1. Basic Constants and Angle Conversions
    c = 299792458;          % Speed of light (m/s)
    lambda = c / f;         % Wavelength (m)
    k = 2 * pi / lambda;    % Free space wave number 
    
    % Convert angles from degrees to radians
    theta1 = deg2rad(theta1_deg);
    theta2 = deg2rad(theta2_deg);
    theta3 = deg2rad(theta3_deg);
    
    % --- Model Validity Check ---
    % The model is not accurate for low grazing angles.
    if theta1 >= pi/2 || theta2 >= pi/2
        warning('Kirchhoff model is inaccurate for grazing angles (theta1 or theta2 >= 90 degrees).');
        sigma_0 = 0;
        return;
    end
    if cos(theta1) == 0 || (cos(theta1) + cos(theta2)) == 0
        warning('Division by zero in geometric factor F. Angles are invalid for this model.');
        sigma_0 = 0;
        return;
    end

    %% 2. Calculate Fresnel Reflection Coefficient (Eq. 9)
    
    % The paper uses theta1 for the Fresnel calculation.
    cos_theta1 = cos(theta1);
    
    if strcmpi(polarization, 'TE')
        % n_eff for TE polarization 
        n_eff = n_complex * cos_theta1;
    elseif strcmpi(polarization, 'TM')
        % n_eff for TM polarization 
        n_eff = n_complex / cos_theta1;
    else
        error('Polarization must be ''TE'' or ''TM''.');
    end
    
    % Fresnel reflection coefficient r [cite: 89]
    r = (1 - n_eff) / (1 + n_eff);
    
    % Power reflection coefficient <rr*> [cite: 93]
    R_fresnel = abs(r)^2;

    %% 3. Calculate Kirchhoff Model Components (Eq. 3, 5-8)

    % Scattering vector components vx, vy (Eq. 6, 7)
    vx = k * (sin(theta1) - sin(theta2) * cos(theta3)); % [cite: 79]
    vy = k * (-sin(theta2) * sin(theta3)); % [cite: 80]
    
    % Rayleigh roughness factor g (Eq. 8)
    g = (k * sigma * (cos(theta1) + cos(theta2)))^2; % [cite: 81]
    
    % Geometric factor F (Eq. 3)
    cos_theta2 = cos(theta2);
    F = (1 + cos_theta1 * cos_theta2 - sin(theta1) * sin(theta2) * cos(theta3)) ...
        / (cos_theta1 * (cos_theta1 + cos_theta2)); % [cite: 65]
    
    % Specular term rho_0 (Eq. 5)
    % The paper uses sinc(x). We assume this is the unnormalized
    % sinc function, sinc(x) = sin(x)/x.
    % MATLAB's sinc is sinc(x) = sin(pi*x)/(pi*x).
    % So, sinc_paper(x) = sinc(x/pi) in MATLAB.
    rho_0 = sinc(vx * lx / pi) * sinc(vy * ly / pi); % [cite: 76]
    rho_0_squared = rho_0^2;
    
    %% 4. Calculate the Diffuse Summation (from Eq. 4)
    
    diffuse_sum = 0;
    M = 50;
    % M = 50; % Truncation limit for the infinite series.
    %         % Increase if g is very large.
    
    v_mag_sq = vx^2 + vy^2; % (vx^2 + vy^2)
    
    % Avoid numerical overflow/underflow by computing the sum
    % iteratively and stopping when terms become negligible.
    g_m = 1;      % g^m
    m_fact = 1;   % m!
    
    for m = 1:M
        g_m = g_m * g;         % g^m
        m_fact = m_fact * m;   % m! (Note: factorial(m) is equivalent)
        
        % This is the m-th term in the series
        exp_term = exp(-(v_mag_sq * L^2) / (4 * m));
        term_m = (g_m / (m_fact * m)) * exp_term; % 
        
        diffuse_sum = diffuse_sum + term_m;
        
        % Convergence check (optional, but good practice)
        if term_m < 1e-12 * diffuse_sum && m > 5
            break;
        end
    end
    
    %% 5. Calculate <pp*>_inf for an infinitely conductive surface (Eq. 4)
    A = lx * ly; % Illuminated area [cite: 63]
    
    diffuse_part = (pi * L^2 * F^2 / A) * diffuse_sum;
    
    % Combine specular and diffuse parts, attenuated by exp(-g)
    rho_inf = exp(-g) * (rho_0_squared + diffuse_part); % 
    
    %% 6. Calculate final <pp*> for finitely conducting surface (Eq. 10)
    
    % Scale the "perfect conductor" result by the Fresnel power coefficient
    sigma_0 = R_fresnel * rho_inf; % [cite: 93]
    
end
% MATLAB's built-in sinc(x) is sin(pi*x)/(pi*x).
% The paper's (and standard physics) sinc(x) is sin(x)/x.
% We use sinc(x/pi) to achieve the paper's definition.