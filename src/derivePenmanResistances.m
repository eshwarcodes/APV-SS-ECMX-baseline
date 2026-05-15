function [r_s_pm, r_a_pm] = derivePenmanResistances(gsCO2, ra_photo, rb_photo, Ts, Pre)
% Converts photosynthesis-level conductance and resistances to Penman–Monteith style
%
% Inputs:
%   gsCO2     - stomatal conductance to CO₂ [µmol CO₂ m⁻² s⁻¹]
%   ra_photo  - aerodynamic resistance from photosynthesis model [m²·s/µmol]
%   rb_photo  - boundary layer resistance from photosynthesis model [m²·s/µmol]
%   Ts        - leaf (or surface) temperature [°C]
%   Pre       - atmospheric pressure [Pa]
%
% Outputs:
%   r_s_pm    - surface resistance [s/m]
%   r_a_pm    - aerodynamic resistance [s/m]

    % Convert gsCO2 [µmol/m²/s] to gsH2O [µmol/m²/s]
    gsH2O = 1.6 * gsCO2;

    % Avoid divide-by-zero
    gsH2O(gsH2O <= 1e-6) = 1e-6;

    % Compute r_s (surface resistance) in [s/m]
    r_s_pm = 1e6 ./ gsH2O;

    % Convert ra + rb to [s/m] using ideal gas conversion
    Pre0 = 101325;     % Standard pressure [Pa]
    Tf = 273.15;       % Standard temp offset
    conversion = (0.0224 * (Ts + 273.15) * Pre0) / (Tf * Pre);  % [µmol/m³ → mol/m³]
    
    % Compute r_a (aerodynamic resistance) in [s/m]
    r_a_pm = (ra_photo + rb_photo) * conversion;
end
