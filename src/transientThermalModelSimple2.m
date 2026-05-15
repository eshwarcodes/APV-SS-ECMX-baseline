function [t_sol, x_sol] = transientThermalModelSimple2(tspan, numGroundX, numGroundY, T_air, RH, wind_speed, LAI_t, gamma, lambda, C_pv_front, C_pv_back, C_ground_top, k_ground, ground_albedo, C_crop, k_crop, C_ground_inner, h_conv_inner, wavelengthDir, L_PV, W_PV, Pr, numPanels, sigma, epsilon_pv, Wg, Front_Irradiance_abs, Rear_Irradiance_abs, solar_rad, T_sky_K_day, T_PV, T_ground, T_crop, VF_accum_front, VF_accum_rear, length_day, f_crop, ri, qy, ab, Ca, Csl, ra, rb, Pre, CT, T0, Vmax, Oa, g1, go, rjv, theta, alpha, cf, c_lar)

N = numGroundX * numGroundY;
ri1 = ri(:,1);
ri9 = ri(:,9);
irradianceCrop_vec = solar_rad(:) * 1e6 / 86400;
R_s_down = solar_rad * 1e6 / 86400;

idx_pv_front = 1:numPanels;
idx_pv_back = numPanels+1:2*numPanels;
idx_ground = 2*numPanels+1:2*numPanels+N;
idx_crop = idx_ground(end)+1:idx_ground(end)+N;
idx_ground_inner = idx_crop(end)+1:idx_crop(end)+2;
total_states = idx_ground_inner(2);
x0 = T_air * ones(total_states,1);

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, ...
                 'OutputFcn', @stopAfterMinutes);

J = speye(total_states);  % Identity mass matrix for future use

[t_sol, x_sol] = ode15s(@thermalODE, tspan, x0, options);

function dxdt = thermalODE(t,x)
    Wg_vec = Wg(:);
    dxdt = zeros(size(x));
    T_amb = T_air + 273.15;
    T_amb_C = T_air;

    % Precompute row-col maps
    persistent row_map col_map
    if isempty(row_map)
        [row_map, col_map] = ind2sub([numGroundY, numGroundX], (1:N)');
    end

    % Extract states
    T_ground_vec = x(idx_ground);
    T_crop_vec   = x(idx_crop);
    T_crop_grid  = reshape(T_crop_vec, numGroundY, numGroundX);

    % --- PV PANEL HEAT BALANCE ---
    T_pv_front_K = x(idx_pv_front) + 273.15;
    T_pv_back_K  = x(idx_pv_back)  + 273.15;
    h_conv_pv_front = estimateConvCoeff(wind_speed, L_PV, W_PV, T_pv_front_K - 273.15, Pr, numPanels);
    h_conv_pv_back  = estimateConvCoeff(wind_speed, L_PV, W_PV, T_pv_back_K - 273.15, Pr, numPanels);

    T_pv_front_K4 = T_pv_front_K.^4;
    T_pv_back_K4  = T_pv_back_K.^4;
    T_g_vec_K4 = (1 - f_crop) .* (T_ground_vec + 273.15).^4 + f_crop .* (T_crop_vec + 273.15).^4;

    VF_front_sum = sum(VF_accum_front, 1)';
    VF_rear_sum  = sum(VF_accum_rear, 1)';

    Q_rad_front = epsilon_pv * sigma * ((VF_accum_front' * T_g_vec_K4) - VF_front_sum .* T_pv_front_K4);
    Q_rad_rear  = epsilon_pv * sigma * ((VF_accum_rear'  * T_g_vec_K4) - VF_rear_sum  .* T_pv_back_K4);

    dxdt(idx_pv_front) = (1 ./ C_pv_front) .* (Front_Irradiance_abs' - h_conv_pv_front .* (x(idx_pv_front) - T_amb_C) - Q_rad_front/24);
    dxdt(idx_pv_back)  = (1 ./ C_pv_back)  .* (Rear_Irradiance_abs'  - h_conv_pv_back  .* (x(idx_pv_back)  - T_amb_C) - Q_rad_rear/24);

    % --- GROUND SURFACE (SOIL) ---
    f_light = exp(-0.7 * c_lar * Wg_vec);
    I_local = irradianceCrop_vec(:) .* f_light;
    R_net_inten = (1 - ground_albedo) .* I_local;
    h_conv_ground = estimateConvCoeff(wind_speed, 0.5, 0.25, T_ground_vec, Pr, N);

    T_crop_lin = T_crop_grid(sub2ind([numGroundY, numGroundX], row_map, col_map));
    T_soil_lin = T_ground_vec;

    rad_contrib = (VF_accum_front * T_pv_front_K4 - sum(VF_accum_front, 2) .* T_g_vec_K4) * epsilon_pv * sigma + ...
                  (VF_accum_rear  * T_pv_back_K4  - sum(VF_accum_rear,  2) .* T_g_vec_K4) * epsilon_pv * sigma;

    dxdt(idx_ground) = (1 ./ C_ground_top) .* ( ...
        R_net_inten - h_conv_ground .* (T_soil_lin - T_amb_C) - k_ground .* (T_soil_lin - x(idx_ground_inner(end))) + ...
        k_crop .* (T_crop_lin - T_soil_lin) + rad_contrib);

    % --- CROP CANOPY ---
    h_conv_crop = estimateConvCoeff(wind_speed, 0.3, 0.1, T_crop_vec, Pr, N);
    row_col_idx = (col_map - 1) * numGroundY + row_map;
    T_crop_px_vec = T_crop_grid(row_col_idx);

    rs = 25;  % surface resistance
    R_s_down_lin = R_s_down(row_col_idx);
    R_n = computeNetRadiation(R_s_down_lin, T_amb_C, RH, T_crop_px_vec, ground_albedo, cf,LAI_t);
    Q_ET = estimateEvapotranspiration(T_crop_px_vec + 273.15, T_amb, RH, wind_speed, LAI_t(row_col_idx), gamma, lambda, R_n, rs, ra);

    T_soil_px_vec = T_ground_vec;
    dxdt(idx_crop) = (1 ./ C_crop) .* ( ...
        R_n - h_conv_crop .* (T_crop_px_vec - T_amb_C) - k_crop .* (T_crop_px_vec - T_soil_px_vec) - Q_ET);

    % --- INNER GROUND ---
    T_ground_avg = mean(T_ground_vec);
    dxdt(idx_ground_inner) = (1 / C_ground_inner) * ( ...
        k_ground * (T_ground_avg - x(idx_ground_inner)) - h_conv_inner * (x(idx_ground_inner) - T_amb_C));
end

end
