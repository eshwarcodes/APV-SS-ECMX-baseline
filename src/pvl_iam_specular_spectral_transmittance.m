function [Tloss_Beam, Tloss_Iso, Tloss_Albedo] = pvl_iam_specular_spectral_transmittance(SurfTilt, AOI, T_stack_allangles)
% Inputs:
% SurfTilt          - scalar or vector tilt in degrees (0–180)
% AOI               - scalar or vector AOI in degrees (0–89)
% T_stack_allangles - 90×341 transmittance matrix [AOI 0°–89° x 300–2000nm in 5nm steps]

% Outputs:
% Tloss_Beam    - 341 x N array of spectral transmittance at specified AOI(s)
% Tloss_Iso     - 341 x 1 spectral transmittance for isotropic diffuse incidence
% Tloss_Albedo  - 341 x 1 spectral transmittance for albedo incidence

% Ensure column vectors
SurfTilt = SurfTilt(:);
AOI = AOI(:);

% Check input dimensions
assert(size(T_stack_allangles, 1) == 90, 'Expected 90 AOI rows in T_stack_allangles.');
assert(size(T_stack_allangles, 2) == 341, 'Expected 341 wavelength columns (300–2000nm).');

% --- 1. Beam transmittance: AOI-specific
n = max(numel(AOI), numel(SurfTilt));  % Match lengths
Tloss_Beam = zeros(341, n);
for i = 1:n
    this_AOI = round(AOI(min(i, end)));   % Allow scalar broadcast

    if this_AOI >= 90
        Tloss_Beam(:, i) = 0;  % Force zero transmittance at 90°
    else
        this_AOI = max(0, min(89, this_AOI)); % Clamp to [0, 89]
        Tloss_Beam(:, i) = T_stack_allangles(this_AOI + 1, :)';  % 341x1
    end
end

% --- 2. Diffuse isotropic: Cosine-weighted average over AOI 0–89°
theta = (0:89)';
theta_rad = deg2rad(theta);
iso_weights = cos(theta_rad) .* sin(theta_rad);
iso_weights = iso_weights / sum(iso_weights);  % Normalize

Tloss_Iso = T_stack_allangles' * iso_weights;  % 341x1

% --- 3. Albedo: Sine-weighted average over 0–89°
albedo_weights = sin(2 * theta_rad);
albedo_weights = albedo_weights / sum(albedo_weights);

Tloss_Albedo = T_stack_allangles' * albedo_weights;  % 341x1

end
