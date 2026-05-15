function R_n = computeNetRadiation(R_s_down, T_air, RH, T_crop, albedo, cf, LAI)
% computeNetRadiation - Net radiation at crop surface with LAI influence
%
% Inputs:
%   R_s_down - [n×1] Incoming shortwave radiation [W/m²]
%   T_air    - [n×1] Air temperature [°C]
%   RH       - [n×1] Relative humidity [%]
%   T_crop   - [n×1] Crop surface temperature [°C]
%   albedo   - [n×1 or scalar] Surface reflectivity
%   cf       - [n×1 or scalar] Cloud fraction
%   LAI      - [n×1 or scalar] Leaf area index
%
% Output:
%   R_n      - [n×1] Net radiation [W/m²]

    % Constants
    sigma = 5.67e-8;  % Stefan-Boltzmann [W/m²/K⁴]
    k_ext = 0.7;      % Light extinction coefficient for fAPAR

    % LAI-based light absorption
    fAPAR = 1 - exp(-k_ext * LAI);        % Fraction of absorbed PAR
    fAPAR = min(max(fAPAR, 0), 1);        % Clamp to [0,1]

    % 1. Wet-bulb temperature approx
    Twb = (T_air .* atan(0.151977 * sqrt(RH + 8.313659))) + ...
          atan(T_air + RH) - atan(RH - 1.676331) + ...
          (0.00391838 * RH.^1.5 .* atan(0.023101 * RH)) - 4.686035;

    % 2. Clear-sky emissivity
    eclear = 0.787 + 0.7641 * log((Twb + 273) ./ 273);

    % 3. Effective sky emissivity
    esky = (1 + 0.0224 * cf - 0.0035 * cf.^2 + 0.00028 * cf.^3) .* eclear;

    % 4. Sky temp [K]
    T_sky_K = (T_air + 273.5) .* esky.^0.25;

    % 5. Crop temp [K]
    T_crop_K = T_crop + 273.15;

    % 6. Shortwave
    R_s_abs = R_s_down;           % absorbed shortwave
    R_s_up  = albedo .* R_s_down;          % reflected portion (not fAPAR-scaled)
    
    % 7. Longwave
    e_s = 0.6108 .* exp((17.27 .* T_air) ./ (T_air + 237.3));
    e_a = (RH ./ 100) .* e_s;
    epsilon_sky = 0.34 - 0.14 * sqrt(e_a);
    
    R_l_down = epsilon_sky .* sigma .* T_sky_K.^4;  % full sky down
    R_l_up   = sigma .* T_crop_K.^4;                % full canopy up


    % 8. Net radiation
    R_n = ((R_s_abs.*fAPAR) - R_s_up) + (epsilon_sky .* sigma .* ((T_sky_K.^4)-(T_crop_K.^4)));

end
