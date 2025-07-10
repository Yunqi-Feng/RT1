function [utd_loss_dB] = calculate_utd_loss(d1, d2, D)
%CALCULATE_UTD_LOSS Calculates the diffraction loss in dB from the UTD
%diffraction coefficient.
%
% Inputs:
%   d1 - Distance from the source to the diffraction point.
%   d2 - Distance from the diffraction point to the observation point.
%   D  - The appropriate UTD diffraction coefficient (D_s or D_h).
%
% Output:
%   utd_loss_dB - The diffraction loss in dB.

    % The spreading factor A(d)
    A_d = sqrt(d1 / (d2 * (d1 + d2)));
    
    % The diffracted field is E_d = E_i * D * A_d
    % The loss is the ratio of the diffracted field to the incident field
    loss_factor = abs(D * A_d);
    
    % Convert to dB
    utd_loss_dB = -20 * log10(loss_factor);

end