function results = run_uncertainty(N, sim_cfg, perturb_base, save_dir, baseline)
% RUN_UNCERTAINTY  Monte-Carlo uncertainty analysis for APV crop simulation.
%
%   results = run_uncertainty(N, sim_cfg, perturb_base, save_dir, baseline)
%
%   N            - number of Monte-Carlo iterations
%   sim_cfg      - struct of constant parameters for simulateCoupledCropThermal4
%   perturb_base - struct of base values for parameters that are perturbed
%   save_dir     - directory to save results (timestamped .mat file)
%   baseline     - (optional) struct with fields for downstream metrics:
%                    .E_APVSS_kWh  — baseline PV energy [kWh]
%                    .CV_ePAR      — baseline spatial CV of ePAR [%]
%                    .CV_IR        — baseline spatial CV of irradiance [%]
%
%   Supports both lettuce and tomato via sim_cfg.crop.type.
%
%   Per-iteration PV energy is reconstructed as:
%     E_i = baseline.E_APVSS_kWh * scales.rad(i)
%   since PV output scales linearly with radiative input.

if nargin < 5, baseline = struct(); end

crop_type = lower(sim_cfg.crop.type);
fprintf('\nRunning %s uncertainty, N = %d\n', crop_type, N);

%% ========== Generate scale factors ==========
rng(1, 'combRecursive');

% Photosynthetic (+-10%)
s_Vcmax  = 1 + 0.10*(2*rand(N,1)-1);
s_alpha  = 1 + 0.10*(2*rand(N,1)-1);
s_g1     = 1 + 0.10*(2*rand(N,1)-1);
s_rjv    = 1 + 0.10*(2*rand(N,1)-1);
s_thetaF = 1 + 0.05*(2*rand(N,1)-1);

% Thermal (+-20%)
s_conv = 1 + 0.20*(2*rand(N,1)-1);

% Structural (+-10%, +-5%)
s_c_lar   = 1 + 0.10*(2*rand(N,1)-1);
s_c_gr    = 1 + 0.10*(2*rand(N,1)-1);
s_theta_m = 1 + 0.05*(2*rand(N,1)-1);
s_DMF     = 1 + 0.05*(2*rand(N,1)-1);
s_SLAmin  = 1 + 0.10*(2*rand(N,1)-1);
s_SLAmax  = 1 + 0.10*(2*rand(N,1)-1);
s_SLAexp  = 1 + 0.10*(2*rand(N,1)-1);
s_betaT   = 1 + 0.10*(2*rand(N,1)-1);
s_betaC   = 1 + 0.10*(2*rand(N,1)-1);

% Canopy optics (+-10%)
s_kext = 1 + 0.10*(2*rand(N,1)-1);

% Radiative input (+-5%)
s_rad = 1 + 0.05*(2*rand(N,1)-1);

% Store scales for reproducibility
scales = struct( ...
    'Vcmax',s_Vcmax, 'alpha',s_alpha, 'g1',s_g1, 'rjv',s_rjv, ...
    'thetaF',s_thetaF, 'conv',s_conv, ...
    'c_lar',s_c_lar, 'c_gr',s_c_gr, 'theta_m',s_theta_m, 'DMF',s_DMF, ...
    'SLAmin',s_SLAmin, 'SLAmax',s_SLAmax, 'SLAexp',s_SLAexp, ...
    'betaT',s_betaT, 'betaC',s_betaC, ...
    'kext',s_kext, 'rad',s_rad);

%% ========== Pre-allocate results ==========
summary_fields = {'meanYield','stdYield','minYield','maxYield','CVYield', ...
                  'harvestMean','harvestStd','harvestIQR', ...
                  'WUE','dTcan','tempMean','tempStd','tempCV'};
empty_summary = cell2struct(repmat({NaN}, numel(summary_fields), 1), ...
                            summary_fields, 1);
results = repmat(empty_summary, N, 1);

%% ========== Unpack sim_cfg for parfor broadcast ==========
sc_tspan       = sim_cfg.tspan;
sc_numGroundX  = sim_cfg.numGroundX;
sc_numGroundY  = sim_cfg.numGroundY;
sc_Tao_daily   = sim_cfg.Tao_daily;
sc_rhao_daily  = sim_cfg.rhao_daily;
sc_wind_daily  = sim_cfg.wind_daily;
sc_gamma       = sim_cfg.gamma;
sc_lambda      = sim_cfg.lambda;
sc_C_pv_front  = sim_cfg.C_pv_front;
sc_C_pv_back   = sim_cfg.C_pv_back;
sc_C_ground_top    = sim_cfg.C_ground_top;
sc_k_ground        = sim_cfg.k_ground;
sc_ground_albedo   = sim_cfg.ground_albedo_int;
sc_C_crop          = sim_cfg.C_crop;
sc_k_crop          = sim_cfg.k_crop;
sc_C_ground_inner  = sim_cfg.C_ground_inner;
sc_wavelengthDir   = sim_cfg.wavelengthDir;
sc_panelWidth      = sim_cfg.panelWidth;
sc_panelHeight     = sim_cfg.panelHeight;
sc_Pr              = sim_cfg.Pr;
sc_numPanels       = sim_cfg.numPanels;
sc_sig             = sim_cfg.sig;
sc_epsilon_pv      = sim_cfg.epsilon_pv;
sc_ri              = sim_cfg.ri;
sc_sp              = sim_cfg.sp;
sc_target_head_weight = sim_cfg.target_head_weight;
sc_planting_density   = sim_cfg.planting_density;
sc_usable_fraction    = sim_cfg.usable_fraction;
sc_max_days           = sim_cfg.max_days;
sc_Ca   = sim_cfg.Ca;
sc_Csl  = sim_cfg.Csl;
sc_Oa   = sim_cfg.Oa;
sc_ra   = sim_cfg.ra;
sc_rb   = sim_cfg.rb;
sc_Pre  = sim_cfg.Pre;
sc_CT   = sim_cfg.CT;
sc_T0   = sim_cfg.T0;
sc_g0   = sim_cfg.g0;
sc_Ws_init  = sim_cfg.Ws_init;
sc_Wg_init  = sim_cfg.Wg_init;
sc_GDD_init = sim_cfg.GDD_init;
sc_c_shoot  = sim_cfg.c_shoot;
sc_c_root   = sim_cfg.c_root;
sc_f_root   = sim_cfg.f_root;
sc_cp       = sim_cfg.cp;
sc_Tbase    = sim_cfg.Tbase;
sc_T_sky_K_day          = sim_cfg.T_sky_K_day;
sc_hours_per_day        = sim_cfg.hours_per_day;
sc_hourNumbersSorted_VF = sim_cfg.hourNumbersSorted_VF;
sc_VF_front  = sim_cfg.viewFactor_GroundToPVFront;
sc_VF_rear   = sim_cfg.viewFactor_GroundToPVRear;
sc_cf_matrix = sim_cfg.cf_matrix;
sc_crop      = sim_cfg.crop;

%% ========== Unpack perturb_base ==========
pb_Vcmax0    = perturb_base.Vcmax0;
pb_alpha     = perturb_base.alpha;
pb_g1        = perturb_base.g1;
pb_rjv       = perturb_base.rjv;
pb_theta     = perturb_base.theta;
pb_h_conv    = perturb_base.h_conv_inner;
pb_c_lar     = perturb_base.c_lar;
pb_c_gr_max  = perturb_base.c_gr_max;
pb_theta_m   = perturb_base.theta_m;
pb_DMF       = perturb_base.dry_matter_fraction;
pb_SLAmin    = perturb_base.SLAmin;
pb_SLAmax    = perturb_base.SLAmax;
pb_SLAexp    = perturb_base.SLAexp;
pb_betaT     = perturb_base.betaT;
pb_betaC     = perturb_base.betaC;
pb_k_ext     = perturb_base.k_ext;
pb_ePAR      = perturb_base.daily_ePAR_totals;
pb_totals    = perturb_base.daily_totals;
pb_moles     = perturb_base.daily_ePAR_mole_totals;
pb_front_abs = perturb_base.daily_front_abs_totals;
pb_rear_abs  = perturb_base.daily_rear_abs_totals;

%% ========== parfor ==========
parfor i = 1:N

    % --- Apply perturbations (scalars) ---
    p_Vcmax0   = pb_Vcmax0   * s_Vcmax(i);
    p_alpha    = pb_alpha    * s_alpha(i);
    p_g1       = pb_g1       * s_g1(i);
    p_rjv      = pb_rjv      * s_rjv(i);
    p_theta    = pb_theta    * s_thetaF(i);
    p_h_conv   = pb_h_conv   * s_conv(i);
    p_c_lar    = pb_c_lar    * s_c_lar(i);
    p_c_gr_max = pb_c_gr_max * s_c_gr(i);
    p_theta_m  = pb_theta_m  * s_theta_m(i);
    p_DMF      = pb_DMF      * s_DMF(i);
    p_SLAmin   = pb_SLAmin   * s_SLAmin(i);
    p_SLAmax   = pb_SLAmax   * s_SLAmax(i);
    p_SLAexp   = pb_SLAexp   * s_SLAexp(i);
    p_betaT    = pb_betaT    * s_betaT(i);
    p_betaC    = pb_betaC    * s_betaC(i);
    p_k_ext    = pb_k_ext    * s_kext(i);

    % --- Apply perturbations (cell arrays — radiative envelope) ---
    rad_i = s_rad(i);
    p_ePAR  = cellfun(@(x) x * rad_i, pb_ePAR,  'UniformOutput', false);
    p_tots  = cellfun(@(x) x * rad_i, pb_totals, 'UniformOutput', false);
    p_moles = cellfun(@(x) x * rad_i, pb_moles,  'UniformOutput', false);
    p_front = pb_front_abs * rad_i;
    p_rear  = pb_rear_abs  * rad_i;

    % --- Build crop struct for this iteration ---
    crop_i = sc_crop;
    if strcmp(crop_type, 'tomato') && isfield(crop_i, 'params')
        crop_i.params.betaT  = p_betaT;
        crop_i.params.betaC  = p_betaC;
        crop_i.params.SLAmin = p_SLAmin;
        crop_i.params.SLAmax = p_SLAmax;
    end

    % --- Call simulation ---
    [~, md] = simulateCoupledCropThermal4( ...
        sc_tspan, sc_numGroundX, sc_numGroundY, ...
        sc_Tao_daily, sc_rhao_daily, sc_wind_daily, ...
        sc_gamma, sc_lambda, ...
        sc_C_pv_front, sc_C_pv_back, ...
        p_front, p_rear, ...
        sc_C_ground_top, sc_k_ground, sc_ground_albedo, sc_C_crop, ...
        sc_k_crop, sc_C_ground_inner, p_h_conv, ...
        sc_wavelengthDir, sc_panelWidth, sc_panelHeight, sc_Pr, ...
        sc_numPanels, sc_sig, sc_epsilon_pv, ...
        sc_ri, sc_sp, p_ePAR, ...
        sc_target_head_weight, sc_planting_density, ...
        p_DMF, sc_usable_fraction, ...
        sc_max_days, p_theta_m, ...
        sc_Ca, sc_Csl, sc_Oa, sc_ra, sc_rb, sc_Pre, ...
        sc_CT, sc_T0, p_Vcmax0, p_g1, sc_g0, p_rjv, p_theta, p_alpha, ...
        sc_Ws_init, sc_Wg_init, sc_GDD_init, ...
        sc_c_shoot, sc_c_root, sc_f_root, sc_cp, ...
        p_c_gr_max, sc_Tbase, p_c_lar, ...
        sc_T_sky_K_day, sc_hours_per_day, ...
        sc_hourNumbersSorted_VF, ...
        sc_VF_front, ...
        sc_VF_rear, ...
        sc_cf_matrix, crop_i, ...
        p_tots, p_moles, p_k_ext, ...
        p_SLAmin, p_SLAmax, ...
        p_SLAexp, p_betaT, p_betaC);

    % --- Extract scalar summary only (no grid storage) ---
    results(i) = extract_summary_single(md, sc_Tao_daily, sc_planting_density);
end

%% ========== Post-process: per-iteration energy & spatial CVs ==========
if isfield(baseline, 'E_APVSS_kWh')
    E_base = baseline.E_APVSS_kWh;
    for i = 1:N
        results(i).E_kWh = E_base * s_rad(i);
    end
else
    for i = 1:N
        results(i).E_kWh = NaN;
    end
end

% CVePAR and CVIR don't change with uniform rad scaling — copy baseline
if isfield(baseline, 'CV_ePAR')
    for i = 1:N, results(i).CVePAR = baseline.CV_ePAR; end
else
    for i = 1:N, results(i).CVePAR = NaN; end
end
if isfield(baseline, 'CV_IR')
    for i = 1:N, results(i).CVIR = baseline.CV_IR; end
else
    for i = 1:N, results(i).CVIR = NaN; end
end

%% ========== Save ==========
timestamp = datestr(now, 'yyyymmdd_HHMMss');
fname = sprintf('uncertainty_%s_N%d_%s.mat', crop_type, N, timestamp);
save_path = fullfile(save_dir, fname);
save(save_path, 'results', 'scales', 'perturb_base', 'baseline', '-v7.3');
fprintf('Saved: %s\n', save_path);

end
