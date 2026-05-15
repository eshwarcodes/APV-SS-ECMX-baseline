function [growth_results, month_data] = simulateCoupledCropThermal4( ...
    tspan, numGroundX, numGroundY, T_air_series, ...
    rhao_daily, wind_daily, gamma, lambda, ...
    C_pv_front, C_pv_back, ...
    Front_Irradiance_abs, Rear_Irradiance_abs, ...
    C_ground_top, k_ground, ground_albedo, C_crop, ...
    k_crop, C_ground_inner, h_conv_inner, ...
    wavelengthDir, L_PV, W_PV, Pr, ...
    numPanels, sigma, epsilon_pv, ...
    ri, sp, daily_ePAR_totals, target_head_weight, ...
    planting_density, dry_matter_fraction, usable_fraction, ...
    max_days, theta_m, Ca, Csl, Oa, ra, rb, Pre, ...
    CT, T0, Vcmax0, g1, g0, rjv, theta, alpha, ...
    Ws_init, Wg_init, GDD_init, c_shoot, c_root, ...
    f_root, cp, c_gr_max, Tbase, c_lar, ...
    T_sky_K_day, hours_per_day, hourNumbersSorted_VF, ...
    viewFactor_GroundToPVFront, viewFactor_GroundToPVRear, ...
    cf_matrix, crop, daily_totals, daily_ePAR_mole_totals,k_ext_new,...
    SLAmin_new,SLAmax_new,...
        SLAexp_new,betaT_new,betaC_new)

%% ===========================
%  CONSTANTS & INDEX MAPS
% ============================
rs = 25;          % surface resistance
qy = [sp(:,1), sp(:,2)];
ab = [sp(:,1), sp(:,3)];

% Sequential day indices for window-ordered VF hours
doy_vec = ceil(hourNumbersSorted_VF / hours_per_day);
[~, ~, day_indices] = unique(doy_vec, 'stable');
numDays = max(day_indices);

%% ===========================
%  INITIALIZE STATE
% ============================
    % ---------------------------------------------------------------
    % Initialize state based on crop type
    % ---------------------------------------------------------------
    switch lower(crop.type)
        
        case 'lettuce'
            Ws = Ws_init * ones(numGroundY, numGroundX);
            Wg = Wg_init * ones(numGroundY, numGroundX);
            GDD = GDD_init * ones(numGroundY, numGroundX);
            
            harvested        = false(numGroundY, numGroundX);
            harvest_day_grid = NaN(numGroundY, numGroundX);
            
            HeadWeight_prev = zeros(numGroundY, numGroundX);

            crop_results    = cell(max_days,1);
            thermal_outputs = cell(max_days,1);
            headweight_outputs = cell(max_days,1);
            crop_data_outputs  = cell(max_days,1);

            Wg_day_output  = cell(max_days,1);
            Ws_day_output  = cell(max_days,1);
            GDD_day_output = cell(max_days,1);
            ET_mmday_output     = cell(max_days,1);
            ET_Rn_mmday_output  = cell(max_days,1);
        
        
        case 'tomato'
            tom = crop.params;
            
            Wprev   = tom.Wprev_init;
            WFprev  = tom.WFprev_init;
            LAIprev = tom.LAIprev_init;
            Nprev   = tom.Nprev_init;
            
            nL = tom.nL;
            Ns = tom.Ns_init;
            Nl = tom.Nl_init;
            Nf = tom.Nf_init;

            AgeL = zeros(numGroundY, numGroundX, nL);
            AgeF = zeros(numGroundY, numGroundX, nL);
            Rc_prev = tom.Rc_prev_init;
            Cstore  = 0;

            harvested        = false(numGroundY, numGroundX);
            harvest_day_grid = NaN(numGroundY, numGroundX);
            
            HeadWeight_prev = zeros(numGroundY, numGroundX);

            crop_results    = cell(max_days,1);
            thermal_outputs = cell(max_days,1);
            headweight_outputs = cell(max_days,1);
            crop_data_outputs  = cell(max_days,1);

            Wg_day_output  = cell(max_days,1);
            Ws_day_output  = cell(max_days,1);
            GDD_day_output = cell(max_days,1);

            ET_mmday_output     = cell(max_days,1);
            ET_Rn_mmday_output  = cell(max_days,1);

        otherwise
            error('Unknown crop.type = %s', crop.type);
    end


    %% ===========================
    % DAILY LOOP
    % ============================
    for d = 1:min(max_days, numDays)

        day_idx = d;  % sequential index into window arrays

        % -----------------------------
        % Extract daily meteorology
        % -----------------------------
        T_air_day   = T_air_series(day_idx);
        RH_air_day  = rhao_daily(day_idx);
        wind_air_day= wind_daily(day_idx);

        solar_rad_full = daily_totals{day_idx};
        daily_ePAR     = daily_ePAR_totals{day_idx};

        Front_abs_day  = Front_Irradiance_abs(day_idx,:);
        Rear_abs_day   = Rear_Irradiance_abs(day_idx,:);

        % -----------------------------
        % Compute representative 3×3 light bins
        % -----------------------------
        light_map = daily_ePAR;
        lm = light_map(:);
        
        mask = isfinite(lm);
        if ~any(mask)
            error('daily_ePAR_totals{%d} contains no finite values.', day_idx);
        end
        
        nBins = 9;
        vmin = min(lm(mask));
        vmax = max(lm(mask));

        if vmin == vmax
            % trivial case: use 9 evenly spaced indices
            N = numel(lm);
            selected_idx = round(linspace(1, N, nBins)).';
        else
            % stratified sampling
            edges   = linspace(vmin, vmax, nBins+1);
            centers = 0.5*(edges(1:end-1)+edges(2:end));
            bin_id  = discretize(lm, edges);
            
            selected_idx = zeros(nBins,1);
            for k = 1:nBins
                idxs = find(bin_id == k);
                if isempty(idxs)
                    % fallback to nearest available
                    [~, rel] = min(abs(lm(mask) - centers(k)));
                    fid = find(mask);
                    selected_idx(k) = fid(rel);
                else
                    [~, rel] = min(abs(lm(idxs) - centers(k)));
                    selected_idx(k) = idxs(rel);
                end
            end
        end
        
        light_3x3 = solar_rad_full(selected_idx);

        % -----------------------------
        % Crop Type: Compute LAI_3x3 etc
        % -----------------------------
        switch lower(crop.type)
            case 'lettuce'
                Wg_3x3  = Wg(selected_idx);
                LAI_3x3 = c_lar .* Wg_3x3;
                f_crop_3x3 = 1 - exp(-k_ext_new * LAI_3x3);
                LAIprev = c_lar .* Wg;   % update wrapper

            case 'tomato'
                LAI_3x3 = LAIprev(selected_idx);
                f_crop_3x3 = 1 - exp(-k_ext_new * LAI_3x3);

                if exist('T_crop_prev','var')
                    T_est_3x3 = T_crop_prev(selected_idx);
                else
                    T_est_3x3 = T_air_day * ones(9,1);
                end

                CO2ppm = crop.params.CO2;
                SLAmin = SLAmin_new;
                SLAmax = SLAmax_new;
                SLAexp=SLAexp_new;

                SLA_3x3 = (SLAmin + (SLAmax-SLAmin).*exp(-SLAexp.*light_3x3)) ./ ...
                         ((1 + betaT_new.*(24 - T_est_3x3)) .* ...
                          (1 + betaC_new.*(CO2ppm - 350)));

                SLA_3x3 = max(SLA_3x3, 1e-12);

                Wg_3x3 = max(LAI_3x3 ./ SLA_3x3, 0);
        end


        % -----------------------------
        % View factor accumulation for today
        % -----------------------------
        hour_idxs = find(day_indices == day_idx);
        Nhr = length(hour_idxs);

        VF_accum_front = zeros(9, numPanels);
        VF_accum_rear  = zeros(9, numPanels);

        for h = hour_idxs(:).'
            VF_F = viewFactor_GroundToPVFront{h};
            VF_R = viewFactor_GroundToPVRear{h};

            for p = 1:numPanels
                vfF = reshape(VF_F{p}, [], 1);
                vfR = reshape(VF_R{p}, [], 1);

                VF_accum_front(:,p) = VF_accum_front(:,p) + vfF(selected_idx);
                VF_accum_rear(:,p)  = VF_accum_rear(:,p)  + vfR(selected_idx);
            end
        end

        if Nhr > 0
            VF_accum_front = VF_accum_front / Nhr;
            VF_accum_rear  = VF_accum_rear  / Nhr;
        end


        % -----------------------------
        % Run 3×3 transient thermal model
        % -----------------------------
        T_PV     = T_air_day * ones(1, numPanels);
        T_ground = T_air_day * ones(1, 9);
        T_crop   = T_air_day * ones(1, 9);

        if exist('Wg_3x3','var')
            Wcrop_input = Wg_3x3;
        else
            Wcrop_input = zeros(9,1);
        end

        [~, x_sol] = transientThermalModelSimple2( ...
            tspan, 3, 3, ...
            T_air_day, RH_air_day, wind_air_day, ...
            LAI_3x3, gamma, lambda, ...
            C_pv_front, C_pv_back, ...
            C_ground_top, k_ground, ground_albedo, ...
            C_crop, k_crop, C_ground_inner, h_conv_inner, ...
            wavelengthDir, L_PV, W_PV, Pr, numPanels, sigma, epsilon_pv, ...
            Wcrop_input, ...
            Front_abs_day, Rear_abs_day, ...
            reshape(light_3x3, [3,3]), ...
            T_sky_K_day(day_idx), ...
            T_PV, T_ground, T_crop, ...
            VF_accum_front, VF_accum_rear, Nhr, ...
            f_crop_3x3, ...
            ri, qy, ab, ...
            Ca, Csl, ra, rb, Pre, ...
            CT, T0, Vcmax0, Oa, ...
            g1, g0, rjv, theta, alpha, ...
            mean(cf_matrix((day_idx-1)*24+1:day_idx*24)), ...
            c_lar);

        T_crop_final = x_sol(end, 28:36);
        T_soil_final = x_sol(end, 19:27);

        coeff_crop = polyfit(light_3x3(:), T_crop_final(:), 1);
        coeff_soil = polyfit(light_3x3(:), T_soil_final(:), 1);

        T_crop_full = coeff_crop(1)*solar_rad_full + coeff_crop(2);
        T_soil_full = coeff_soil(1)*solar_rad_full + coeff_soil(2);


        % -----------------------------
        % Daily Crop Growth Update
        % -----------------------------
        switch lower(crop.type)

            case 'lettuce'
                Rn = computeNetRadiation(solar_rad_full * 1e6/86400, ...
                                         T_air_day, RH_air_day, ...
                                         T_crop_full, ground_albedo, ...
                                         mean(cf_matrix((day_idx-1)*24+1:day_idx*24)), ...
                                         LAIprev);

                ET_mmday = estimateET_mmday(T_crop_full+273.15, T_air_day+273.15, ...
                                            RH_air_day, wind_air_day, ...
                                            LAIprev, gamma, lambda, ...
                                            Rn, rs, ra);
                
                [~, HW_today, Wg_new, Ws_new, GDD_new, crop_data] = ...
                    simulateLettuceGrowthSingleDay( ...
                        ri, sp, T_crop_full, daily_ePAR_totals{day_idx}, ...
                        target_head_weight, planting_density, ...
                        dry_matter_fraction, usable_fraction, ...
                        Ca, Csl, Oa, ra, rb, Pre, RH_air_day, ...
                        CT, T0, Vcmax0, g1, g0, rjv, theta, alpha, ...
                        Ws, Wg, GDD, c_shoot, c_root, ...
                        f_root, cp, c_gr_max, Tbase, c_lar);

                % harvest decisions
                newHarvest = (HW_today >= target_head_weight) | ...
                             (GDD_new > theta_m) | ...
                             (d >= max_days);

                newly = newHarvest & ~harvested;
                harvested(newly) = true;
                harvest_day_grid(newly) = d;

                % freeze harvested pixels
                HW_today(harvested) = HeadWeight_prev(harvested);

                Wg(~harvested)  = Wg_new(~harvested);
                Ws(~harvested)  = Ws_new(~harvested);
                GDD(~harvested) = GDD_new(~harvested);

                HeadWeight_prev = HW_today;

                % store outputs
                crop_results{d}        = T_crop_full;
                thermal_outputs{d}     = T_soil_full;
                headweight_outputs{d}  = HW_today;
                crop_data_outputs{d}   = crop_data;
                Wg_day_output{d}       = Wg;
                Ws_day_output{d}       = Ws;
                GDD_day_output{d}      = GDD;
                ET_mmday_output{d}     = ET_mmday;
                ET_Rn_mmday_output{d}  = ET_mmday;


            case 'tomato'
                tom = crop.params;

                DLI_today  = daily_ePAR_mole_totals{day_idx};
                sunrise_t  = tom.sunrise(day_idx);
                sunset_t   = tom.sunset(day_idx);
                Tday_hours = tom.Thours((day_idx-1)*24+1:day_idx*24);

                Rn = computeNetRadiation(solar_rad_full*1e6/86400, ...
                                         T_air_day, RH_air_day, ...
                                         T_crop_full, ground_albedo, ...
                                         mean(cf_matrix((day_idx-1)*24+1:day_idx*24)), ...
                                         LAIprev);

                ET_mmday = estimateET_mmday(T_crop_full+273.15, ...
                                            T_air_day+273.15, ...
                                            RH_air_day, wind_air_day, ...
                                            LAIprev, gamma, lambda, ...
                                            Rn, rs, ra);

                [Wnow, WFnow, LAInow, Nnow, Ns, Nl, Nf, AgeL, AgeF, Rc_prev, Cstore] = ...
                    simulateTomatoGrowthSingleDay_vectorized( ...
                        T_crop_full, DLI_today, Wprev, WFprev, ...
                        LAIprev, Nprev, Ns, Nl, Nf, ...
                        AgeL, AgeF, Rc_prev, d, tom, ...
                        Tday_hours, tom.Smax, tom.eta_store, ...
                        tom.eta_withd, Cstore);

                % harvest rule
                newHarvest = (WFnow >= tom.WF_target) | (d >= max_days);
                newly = newHarvest & ~harvested;

                harvested(newly) = true;
                harvest_day_grid(newly) = d;

                FullFruit = WFnow;
                FullFruit(harvested) = HeadWeight_prev(harvested);

                HeadWeight_prev = FullFruit;

                % roll forward states
                Wprev   = Wnow;
                WFprev  = WFnow;
                LAIprev = LAInow;
                Nprev   = Nnow;

                % record
                crop_results{d}      = T_crop_full;
                thermal_outputs{d}   = T_soil_full;
                headweight_outputs{d}= FullFruit;
                crop_data_outputs{d} = struct('N',Nnow,'Rc',Rc_prev);

                Wg_day_output{d}     = Wnow;
                Ws_day_output{d}     = WFnow;
                GDD_day_output{d}    = LAInow;  % placeholder

                ET_mmday_output{d}     = ET_mmday;
                ET_Rn_mmday_output{d}  = ET_mmday;

                T_crop_prev = T_crop_full;
        end

        if all(harvested(:))
            break;
        end
    end

    % =====================================
    % Summary
    % =====================================
    growth_results = [ ...
        mean(harvest_day_grid(:), 'omitnan'), ...
        mean(headweight_outputs{min(d,min(max_days,numDays))}(:), 'omitnan') ...
    ];

    month_data = struct( ...
        'Wg', {Wg_day_output}, ...
        'Ws', {Ws_day_output}, ...
        'HeadWeight', {headweight_outputs}, ...
        'GDD', {GDD_day_output}, ...
        'harvest_day_grid', harvest_day_grid, ...
        'crop_temperature', {crop_results}, ...
        'soil_temperature', {thermal_outputs}, ...
        'crop_data', {crop_data_outputs}, ...
        'ET_mm_day', {ET_mmday_output}, ...
        'ET_mm_day_Rn', {ET_Rn_mmday_output} ...
    );

end


