function ET = estimateET_mmday( ...
    T_crop_K, T_amb_K, RH, wind_speed, LAI_t, ...
    gamma, lambda, R_n, r_s, r_a)

    % Constants
    k_ext = 0.7;                         % extinction coefficient
    fAPAR = 1 - exp(-k_ext .* LAI_t);    % canopy absorption fraction
    fAPAR = max(min(fAPAR, 1), 0);       % clamp to [0,1]

    % Convert temperatures to Celsius
    T_crop = T_crop_K - 273.15;
    T_amb  = T_amb_K  - 273.15;

    % Saturation vapor pressure at crop temp (kPa)
    e_s = 0.6108 .* exp((17.27 .* T_crop) ./ (T_crop + 237.3));

    % Actual vapor pressure from ambient (kPa)
    e_a = (RH ./ 100) .* 0.6108 .* exp((17.27 .* T_amb) ./ (T_amb + 237.3));

    % Vapor Pressure Deficit (VPD, kPa)
    VPD = max(0, e_s - e_a);

    % Slope of saturation vapor pressure curve at crop temp (kPa/°C)
    delta = 4098 .* e_s ./ (T_crop + 237.3).^2;

    % Handle invalid resistances
    bad_rs = isnan(r_s) | r_s > 1000;
    bad_ra = isnan(r_a) | r_a > 1000;
    if any(bad_rs(:))
        r_s(bad_rs) = max(50 ./ LAI_t(bad_rs), 200);
    end
    if any(bad_ra(:))
        r_a(bad_ra) = 200;
    end

    % Absorbed canopy net radiation
    R_n_canopy = fAPAR .* R_n;

    %% ----- Penman–Monteith ET (mm/day) -----
    ET_pm = (0.408 .* delta .* R_n_canopy + ...
             gamma .* (900 ./ (T_crop + 273)) .* wind_speed .* VPD) ./ ...
            (delta + gamma .* (1 + r_s ./ r_a));
    ET_pm = max(ET_pm,0);
    %% ----- Convert ET to latent heat flux [W/m²] -----
    Q_ET = lambda .* (ET_pm ./ 86400);   % W/m²
    
    %% ----- Priestley–Taylor-style dynamic cap -----
    AlphaMin = 1.0;
    AlphaMax = 1.3;
    alpha = 1.0 + 0.2 .* tanh((VPD - 1.0) ./ 0.5) ...
                 + 0.1 .* tanh((wind_speed - 2.0) ./ 1.5);
    alpha = min(max(alpha, AlphaMin), AlphaMax);

    EF_PT = delta ./ (delta + gamma);          % Priestley–Taylor EF
    Qcap  = alpha .* EF_PT .* R_n_canopy;      % W/m² cap
    Qcap = max(Qcap,0);
    % Apply cap
    Q_ET = min(Q_ET, Qcap);

    %% ----- Final ET in mm/day -----
    ET = (Q_ET ./ lambda) * 86400;

end
