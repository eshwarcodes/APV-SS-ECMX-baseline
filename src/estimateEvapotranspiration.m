function Q_ET = estimateEvapotranspiration(...
    T_crop_K, T_amb_K, RH, wind_speed, LAI_t, ...
    gamma, lambda, R_n, r_s, r_a)

    % Constants
    k_ext = 0.7;                    % extinction coefficient
    fAPAR = 1 - exp(-k_ext .* LAI_t);  % canopy absorption fraction
    fAPAR = max(min(fAPAR, 1), 0);  % clamp to [0,1]

    % Convert temperatures to Celsius
    T_crop = T_crop_K - 273.15;
    T_amb = T_amb_K - 273.15;

    % Saturation vapor pressure (kPa)
    e_s = 0.6108 .* exp((17.27 .* T_crop) ./ (T_crop + 237.3));

    % Actual vapor pressure (kPa)
    e_a = (RH ./ 100) .* ...
        0.6108 .* exp((17.27 .* T_amb) ./ (T_amb + 237.3));

    % Vapor Pressure Deficit (VPD)
    VPD = max(0, e_s - e_a);

    % Slope of saturation vapor pressure curve (delta)
    delta = 4098 .* e_s ./ (T_crop + 237.3).^2;

    % Handle invalid resistances
    r_s(isnan(r_s) | r_s > 1000) = max(50 ./ LAI_t(isnan(r_s) | r_s > 1000), 200);
    r_a(isnan(r_a) | r_a > 1000) = 200;
    R_n_canopy = fAPAR .* R_n; 
    %% ----- 1. Resistance-limited ET (Penman-Monteith) -----
    ET_pm = (0.408 .* delta .* R_n_canopy + ...
             gamma .* (900 ./ (T_crop + 273)) .* wind_speed .* VPD) ./ ...
            (delta + gamma .* (1 + r_s ./ r_a));
    %ET_pm = max(ET_pm, 0);     % mm/day

    %% ----- 2. Energy-limited ET -----
    %R_n_canopy = fAPAR .* R_n;        % absorbed radiation
    %ET_energy = R_n_canopy ./ lambda; % W/m² / J/kg → kg/m²/s = mm/s
    %ET_energy = ET_energy .* 86400;   % convert to mm/day

    %% ----- 3. Final ET (take the limiting value) -----
    ET = min(max(0, ET_pm),6);

    %% Convert final ET to latent heat loss [W/m²]
    Q_ET = lambda .* (ET ./ 86400);   % kg/m²/s → W/m²
    %Q_ET = min(Q_ET, 0.8 .* R_n);       % cap at 80% of Rn

end
