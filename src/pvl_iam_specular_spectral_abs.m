function [A_Beam, A_Iso, A_Albedo] = pvl_iam_specular_spectral_abs(SurfTilt, AOI, A_stack_allangles)
% Inputs:
% SurfTilt          - scalar or vector tilt in degrees (0–180)
% AOI               - scalar or vector AOI in degrees (0–89)
% R_stack_allangles - 90×341 reflectance matrix [AOI from 0°–89° x 300–2000nm in 5nm steps]

% Ensure column vectors
SurfTilt = SurfTilt(:);
AOI = AOI(:);

% Check input dimensions
assert(size(A_stack_allangles, 1) == 90, 'Expected 90 AOI rows in R_stack_allangles.');
assert(size(A_stack_allangles, 2) == 341, 'Expected 341 wavelength columns (300–2000nm).');

% --- 1. Beam reflectance: AOI-specific reflectance (no spectral weighting)
% Clamp AOI between 0 and 89 and grab the corresponding reflectance rows
n = max(numel(AOI), numel(SurfTilt));  % match lengths
A_Beam = zeros(341, n);
for i = 1:n
    this_AOI = round(AOI(min(i,end)));   % allow scalar broadcast
    if(AOI<=89)
    this_AOI = max(0, min(89, this_AOI)); % clamp to 0–89
    A_Beam(:, i) = A_stack_allangles(this_AOI+1, :)'; % 341x1
    else
        A_Beam(:, i) = 0; % 341x1
    end
end

% --- 2. Diffuse isotropic: Cosine-weighted average over AOI 0–89°
theta = (0:89)';  % degrees
theta_rad = deg2rad(theta);
iso_weights = cos(theta_rad) .* sin(theta_rad);
iso_weights = iso_weights / sum(iso_weights);  % normalize
A_Iso = A_stack_allangles' * iso_weights;  % 341×1

% --- 3. Albedo: Sine-weighted average over reflected hemisphere (mirror 0–89°)
albedo_weights = sin(2 * theta_rad);
albedo_weights = albedo_weights / sum(albedo_weights);
A_Albedo = A_stack_allangles' * albedo_weights;  % 341×1

end
