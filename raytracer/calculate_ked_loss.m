function diffraction_loss_dB = calculate_ked_loss(h, d1, d2, freq)
% calculate_ked_loss: Calculates diffraction loss using the Knife-Edge Diffraction model.
%
% Inputs:
%   h:      Height of the obstruction above the line of sight (in meters).
%   d1:     Distance from the source to the obstruction (in meters).
%   d2:     Distance from the obstruction to the destination (in meters).
%   lambda: Wavelength of the signal (in meters).
%
% Output:
%   diffraction_loss_dB: Diffraction loss in dB.

    % Calculate the Fresnel-Kirchhoff diffraction parameter 'v'
    lambda = 3e8 / freq;
    v = h * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));

    % Calculate the diffraction loss in dB
    if v > -0.78
        diffraction_loss_dB = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
    else
        diffraction_loss_dB = 0;
    end
end