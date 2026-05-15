tic
clear; close all; tic;

% Get main folder path
%mainFolder = 'E:\Yuma_APVSS_dual_axis\022226\baseline_radiative';
addpath('E:\Yuma_APVSS_dual_axis\022226\uncertainty\shared')
addpath('E:\Yuma_APVSS_dual_axis\022226\models')
addpath('E:\Yuma_APVSS_dual_axis\022226\uncertainty\lettuce')
mainFolder ='E:\Yuma_APVSS_dual_axis\022226\baseline_radiative';
% Load base plots
load('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\plots_base.mat');
%load('E:\Yuma_opaque_fixed_tilt\Yuma_opaque_Fixed_tilt_p1_72007\plots_base.mat');
for i=1:8760
% Store day length once per day
    if mod(i, 24) == 0
        daylength(i/24) = S0;
        sunrise(i/24)=ltsr;
    end
end
sunrise = sunrise';
sunset=sunrise+daylength;
% Define base scratch directories
%directory_Yuma_APVSS = ...
%    'E:\Yuma_APVSS_dual_axis\022226\baseline_radiative';
directory_Yuma_APVSS = ...
   'E:\Yuma_APVSS_dual_axis\022226\baseline_radiative';
dataFolder = 'E:\Yuma_APVSS_dual_axis\022226\baseline_radiative';

load(fullfile(dataFolder, 'GI_Yuma_APVSS.mat'));
load(fullfile(dataFolder, 'PI_Yuma_APVSS.mat'));
load(fullfile(dataFolder, 'PP_Yuma_APVSS.mat'));
load(fullfile(dataFolder, 'VF_Yuma_APVSS.mat'));
% Ground Irradiance
dataFolder_Yuma_APVSS_GI = fullfile(directory_Yuma_APVSS, 'groundIrradiance');
filePattern = 'irradianceGround_hour_*.mat';
regexPattern = 'irradianceground_hour_(\d+)\.mat';
[groundIrradiance_Yuma_APVSS,hourNumbersSorted_Yuma_APVSS_grIr, sortIdx] = loadMatFilesFromDir(dataFolder_Yuma_APVSS_GI, filePattern, regexPattern);
fprintf('Loaded groundIrradiance from %s (%d entries)\n', dataFolder_Yuma_APVSS_GI, numel(hourNumbersSorted_Yuma_APVSS_grIr));
% Panel Irradiance
dataFolder_Yuma_APVSS = fullfile(directory_Yuma_APVSS, 'panelIrradiance');
filePattern = 'panelirradiance_hour_*.mat';
regexPattern = 'panelirradiance_hour_(\d+)\.mat';
[panelirradiance_Yuma_APVSS,hourNumbersSorted_Yuma_APVSS_pair, sortIdx] = loadMatFilesFromDir(dataFolder_Yuma_APVSS, filePattern, regexPattern);
fprintf('Loaded panelIrradiance from %s (%d entries)\n', dataFolder_Yuma_APVSS, numel(hourNumbersSorted_Yuma_APVSS_pair));
% PV Performance
dataFolder_Yuma_APVSS = fullfile(directory_Yuma_APVSS, 'PVPerformance');
filePattern = 'PVPanelPerformance_hour_*.mat';
regexPattern = 'pvpanelperformance_hour_(\d+)\.mat';
[PV_Yuma_APVSS, hourNumbersSorted_Yuma_APVSS_PV, sortIdx] = loadMatFilesFromDir(dataFolder_Yuma_APVSS, filePattern, regexPattern);
fprintf('Loaded PVPerformance from %s (%d entries)\n', dataFolder_Yuma_APVSS, numel(hourNumbersSorted_Yuma_APVSS_PV));
dataFolder_Yuma_APVSS = fullfile(directory_Yuma_APVSS, 'viewfactor_groundPV');
filePattern = 'viewFactor_hour_*.mat';
regexPattern = 'viewfactor_hour_(\d+)\.mat';

[VF_Yuma_APVSS, hourNumbersSorted_VF, sortIdx] = ...
    loadMatFilesFromDir(dataFolder_Yuma_APVSS, filePattern, regexPattern);
fprintf('Loaded PVPerformance from %s (%d entries)\n', dataFolder_Yuma_APVSS, numel(hourNumbersSorted_VF));

hours_per_day = 24;
numDays       = 365;
secondsPerHour = 3600;

% Map hours → day index
day_indices = ceil(hourNumbersSorted_Yuma_APVSS_grIr / hours_per_day);   % 4284×1

% Preallocate
daily_ePAR_totals = cell(numDays, 1);
daily_totals      = cell(numDays, 1);

for d = 1:numDays

    mask = (day_indices == d);   % logical vector for hours of this day
    hour_list = find(mask);      % indices into groundIrradiance_Yuma_APVSS

    if isempty(hour_list)
        daily_ePAR_totals{d} = zeros(numGroundY, numGroundX);
        daily_totals{d}      = zeros(numGroundY, numGroundX);
        continue
    end

    % Initialize accumulators
    sum_ePAR  = zeros(numGroundY, numGroundX);
    sum_total = zeros(numGroundY, numGroundX);

    for h = hour_list(:)'   % loop over hours in this day

        ePARh = groundIrradiance_Yuma_APVSS.ePAR_inten_all{h};
        PARh  = groundIrradiance_Yuma_APVSS.PAR_inten_all{h};
        IRh   = groundIrradiance_Yuma_APVSS.IR_inten_all{h};

        if isempty(ePARh)
            continue
        end

        sum_ePAR  = sum_ePAR  + (ePARh * secondsPerHour / 1e6);
        sum_total = sum_total + ((PARh + IRh) * secondsPerHour / 1e6);
    end

    daily_ePAR_totals{d} = sum_ePAR;
    daily_totals{d}      = sum_total;
end

% Daily totals per panel (integrated irradiance)
% Map hours → day index
day_indices = ceil(hourNumbersSorted_Yuma_APVSS_pair / hours_per_day);   % 4284×1

numPanels = 9;

daily_front_totals = zeros(numDays, numPanels);
daily_rear_totals  = zeros(numDays, numPanels);

for h = 1:length(day_indices)
    d = day_indices(h);

    front_vals = panelirradiance_Yuma_APVSS.Front_Irradiance_inten_all{h};   % [1×9]
    rear_vals  = panelirradiance_Yuma_APVSS.Rear_Irradiance_inten_all{h};

    if isempty(front_vals)
        continue
    end

    daily_front_totals(d,:) = daily_front_totals(d,:) + ...
        (front_vals' * secondsPerHour / 1e6);

    daily_rear_totals(d,:) = daily_rear_totals(d,:) + ...
        (rear_vals'  * secondsPerHour / 1e6);
end

daily_front_abs_totals = zeros(numDays, numPanels);
daily_rear_abs_totals  = zeros(numDays, numPanels);

for h = 1:length(day_indices)
    d = day_indices(h);

    front_specs = PV_Yuma_APVSS.Front_Irradiance_abs_all{h};   % cell{1×9}, each 1882×1
    rear_specs  = PV_Yuma_APVSS.Rear_Irradiance_abs_all{h};

    if isempty(front_specs)
        continue
    end

    for p = 1:numPanels
        f = front_specs{p};
        r = rear_specs{p};
        if isempty(f), continue; end

        front_W = trapz(wavelengthDir, f);
        rear_W  = trapz(wavelengthDir, r);

        daily_front_abs_totals(d,p) = daily_front_abs_totals(d,p) + ...
                                      front_W * secondsPerHour / 1e6;
        daily_rear_abs_totals(d,p)  = daily_rear_abs_totals(d,p)  + ...
                                      rear_W * secondsPerHour  / 1e6;
    end
end

Twb = (Tao .* atan(0.151977*sqrt(rhao+8.313659))) + ...
      atan(Tao + rhao) - atan(rhao - 1.676331) + ...
      0.00391838 * (rhao.^1.5).*atan(0.023101*rhao) - 4.686035;

eclear = 0.787 + 0.7641 * log((Twb + 273) / 273);

esky = (1 + 0.0224*cf_matrix - 0.0035*(cf_matrix.^2) + ...
       0.00028*(cf_matrix.^3)) .* eclear;

T_sky_K = (Tao + 273.5).* esky.^0.25;

T_sky_K_day = arrayfun(@(d) mean(T_sky_K((d-1)*24+1:d*24)), 1:365).';

numHours = length(VF_Yuma_APVSS.viewFactor_FrontPVGround_sparse_all);

viewFactor_GroundToPVFront = cell(numHours,1);
viewFactor_GroundToPVRear  = cell(numHours,1);

for d = 1:numHours

    panelFactorsFront = VF_Yuma_APVSS.viewFactor_FrontPVGround_sparse_all{d};
    panelFactorsRear  = VF_Yuma_APVSS.viewFactor_RearPVGround_sparse_all{d};

    [VF_G2PV_front, ~] = convertPVToGroundToGroundToPV(...
        panelFactorsFront, panelWidth, panelHeight, ...
        groundXmin, groundXmax, groundYmin, groundYmax, ...
        numGroundX, numGroundY);

    [VF_G2PV_rear, ~] = convertPVToGroundToGroundToPV(...
        panelFactorsRear, panelWidth, panelHeight, ...
        groundXmin, groundXmax, groundYmin, groundYmax, ...
        numGroundX, numGroundY);

    viewFactor_GroundToPVFront{d} = VF_G2PV_front;
    viewFactor_GroundToPVRear{d}  = VF_G2PV_rear;
end



%% === DAILY EPAR IN MOLES (VECTORIZED + CLEAN) ===
hours_per_day = 24;
numDays = 365;
secondsPerHour = 3600;

% Map each hour → day 1..365
day_indices = ceil(hourNumbersSorted_Yuma_APVSS_grIr/ hours_per_day);

% Preallocate daily cell outputs
daily_ePAR_mole_totals = cell(numDays,1);

for d = 1:numDays

    % find hours belonging to day d
    hour_list = find(day_indices == d);

    if isempty(hour_list)
        daily_ePAR_mole_totals{d} = zeros(numGroundY, numGroundX);
        continue
    end

    sum_ePAR_moles = zeros(numGroundY, numGroundX);

    for h = hour_list(:)'   % iterate only hours of this day

        ePARh = groundIrradiance_Yuma_APVSS.ePAR_inten_mmole_all{h};   % µmol m⁻² s⁻¹ at each grid point

        if isempty(ePARh)
            continue
        end

        % Convert to total moles per m² for the hour:
        % µmol/m²/s → mol/m²/hour = (µmol)*3600 / 1e6
        sum_ePAR_moles = sum_ePAR_moles + (ePARh * secondsPerHour / 1e6);
    end

    daily_ePAR_mole_totals{d} = sum_ePAR_moles;
end

% --- Create a time vector corresponding to each hour over the year ---
time_hours = (1:length(Tao))'/24;

% Reshape and average each variable over 24-hour blocks (rows = days)
Tao_daily   = mean(reshape(Tao,      24, []).', 2);  % Temperature [365 x 1]
wind_daily  = mean(reshape(wind_vel, 24, []).', 2);  % Wind speed
rhao_daily  = mean(reshape(rhao,     24, []).', 2);  % Relative humidity
DNI_daily   = mean(reshape(DNI,      24, []).', 2);  % Direct irradiance
DHI_daily   = mean(reshape(DHI,      24, []).', 2);  % Diffuse irradiance
ground_albedo_int = trapz(wavelengthDir,ground_albedo.*dirHorznIrradNormHourlyMod(:,14))/trapz(wavelengthDir,dirHorznIrradNormHourlyMod(:,14));
sig = 5.67*1e-8;
epsilon_pv = 0.85;
k_ext_new=0.7;
SLAmin_new=tomParams.SLAmin;
SLAmax_new=tomParams.SLAmax;
betaT_new=tomParams.betaT;
betaC_new=tomParams.betaC;
SLAexp_new=0.471;
[growth_results_lettuce, month_data_lettuce] = simulateCoupledCropThermal4( ...
  tspan, numGroundX, numGroundY, Tao_daily, ...
  rhao_daily, wind_daily, gamma, lambda, ...
  C_pv_front, C_pv_back, ...
  daily_front_abs_totals, daily_rear_abs_totals, ...
  C_ground_top, k_ground, ground_albedo_int, C_crop, ...
  k_crop, C_ground_inner, h_conv_inner, ...
  wavelengthDir, panelWidth, panelHeight, Pr, ...
  numPanels, sig, epsilon_pv, ...
  ri, sp, daily_ePAR_totals, target_head_weight, ...
  planting_density, dry_matter_fraction, usable_fraction, ...
  max_days, theta_m, Ca, Csl, Oa, ra, rb, Pre, ...
  CT, T0, Vcmax0, g1, g0, rjv, theta, alpha, ...
  Ws_init, Wg_init, GDD_init, c_shoot, c_root, ...
  f_root, cp, c_gr_max, Tbase, c_lar, ...
  T_sky_K_day, hours_per_day, hourNumbersSorted_VF, ...
  viewFactor_GroundToPVFront, viewFactor_GroundToPVRear, cf_matrix, ...
  cropStructLettuce,daily_totals,daily_ePAR_mole_totals,k_ext_new,...
  SLAmin_new,SLAmax_new,...
        SLAexp_new,betaT_new,betaC_new);


toc
save('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\Lettuce_Yuma_APVSS.mat','growth_results_lettuce','month_data_lettuce','-v7.3');    

BaselineInputs = struct();
BaselineInputs.c_lar = c_lar;
BaselineInputs.c_gr_max = c_gr_max;
BaselineInputs.theta_m = theta_m;
BaselineInputs.dry_matter_fraction = dry_matter_fraction;
BaselineInputs.tspan = tspan;
BaselineInputs.numGroundX = numGroundX;
BaselineInputs.numGroundY = numGroundY;

BaselineInputs.Tao_daily = Tao_daily;
BaselineInputs.rhao_daily = rhao_daily;
BaselineInputs.wind_daily = wind_daily;
BaselineInputs.cf_matrix = cf_matrix;

BaselineInputs.daily_front_abs_totals = daily_front_abs_totals;
BaselineInputs.daily_rear_abs_totals  = daily_rear_abs_totals;
BaselineInputs.daily_ePAR_mole_totals = daily_ePAR_mole_totals;
BaselineInputs.daily_ePAR_totals = daily_ePAR_totals;
BaselineInputs.daily_totals = daily_totals;

BaselineInputs.T_sky_K_day = T_sky_K_day;
BaselineInputs.viewFactor_GroundToPVFront = viewFactor_GroundToPVFront;
BaselineInputs.viewFactor_GroundToPVRear  = viewFactor_GroundToPVRear;

BaselineInputs.cropStruct = cropStructLettuce;

% Crop physiology parameters
BaselineInputs.Vcmax0 = Vcmax0;
BaselineInputs.alpha  = alpha;
BaselineInputs.g1     = g1;
BaselineInputs.g0     = g0;
BaselineInputs.rjv    = rjv;
BaselineInputs.theta  = theta;
BaselineInputs.h_conv_inner = h_conv_inner;

save('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\Baseline_Inputs_Lettuce_APVSS.mat','BaselineInputs','-v7.3');

% [growth_results_tomato, month_data_tomato] = simulateCoupledCropThermal4( ...
%   tspan, numGroundX, numGroundY, Tao_daily, ...
%   rhao_daily, wind_daily, gamma, lambda, ...
%   C_pv_front, C_pv_back, ...
%   daily_front_abs_totals, daily_rear_abs_totals, ...
%   C_ground_top, k_ground, ground_albedo_int, C_crop, ...
%   k_crop, C_ground_inner, h_conv_inner, ...
%   wavelengthDir, panelWidth, panelHeight, Pr, ...
%   numPanels, sig, epsilon_pv, ...
%   ri, sp, daily_ePAR_totals, target_head_weight, ...
%   planting_density, dry_matter_fraction, usable_fraction, ...
%   max_days, theta_m, Ca, Csl, Oa, ra, rb, Pre, ...
%   CT, T0, Vcmax0, g1, g0, rjv, theta, alpha, ...
%   Ws_init, Wg_init, GDD_init, c_shoot, c_root, ...
%   f_root, cp, c_gr_max, Tbase, c_lar, ...
%   T_sky_K_day, hours_per_day, hourNumbersSorted_VF, ...
%   viewFactor_GroundToPVFront, viewFactor_GroundToPVRear, cf_matrix, ...
%   cropTom,daily_totals,daily_ePAR_mole_totals);

% Get main folder path
mainFolder = 'E:\Yuma_opaque_fixed_tilt\Yuma_opaque_Fixed_tilt_p1_72007';

% Load base plots
%load('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\plots_base.mat');
%load('E:\Yuma_APVSS_dual_axis\Yuma_APVSS_Dual_axis_p1_72134\plots_base.mat');
% Define base scratch directories
directory_Yuma_APVSS = ...
    'E:\Yuma_opaque_fixed_tilt\Yuma_opaque_Fixed_tilt_p1_72007';

%load('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\plots_base.mat');
%dataFolder_Yuma_APVSS = fullfile(directory_Yuma_APVSS, 'viewfactor_groundPV');
%filePattern = 'viewFactor_hour_*.mat';
%regexPattern = 'viewfactor_hour_(\d+)\.mat';

%[VF_Yuma_APVSS, hourNumbersSorted_VF, sortIdx] = ...
%    loadMatFilesFromDir(dataFolder_Yuma_APVSS, filePattern, regexPattern);

%---- Load radiative-derived frozen data ----
%load('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\Baseline_Inputs_Lettuce.mat');
%load('E:\Yuma_APVSS_dual_axis\022226\baseline_crop\Baseline_Inputs_Lettuce2.mat')
% --- Unpack struct once ---
tspan   = BaselineInputs.tspan;
numGroundX = BaselineInputs.numGroundX;
numGroundY = BaselineInputs.numGroundY;

Tao_daily  = BaselineInputs.Tao_daily;
rhao_daily = BaselineInputs.rhao_daily;
wind_daily = BaselineInputs.wind_daily;

daily_front_abs_totals = BaselineInputs.daily_front_abs_totals;
daily_rear_abs_totals  = BaselineInputs.daily_rear_abs_totals;

daily_ePAR_totals      = BaselineInputs.daily_ePAR_totals;
daily_ePAR_mole_totals = BaselineInputs.daily_ePAR_mole_totals;
daily_totals           = BaselineInputs.daily_totals;

viewFactor_GroundToPVFront = BaselineInputs.viewFactor_GroundToPVFront;
viewFactor_GroundToPVRear  = BaselineInputs.viewFactor_GroundToPVRear;

T_sky_K_day = BaselineInputs.T_sky_K_day;
cf_matrix   = BaselineInputs.cf_matrix;
% ---- Define constants explicitly ----
gamma_var = 0.0660;
sig = 5.67e-8;
epsilon_pv = 0.85;
hours_per_day = 24;
ground_albedo_int = trapz(wavelengthDir,ground_albedo.*dirHorznIrradNormHourlyMod(:,14))/trapz(wavelengthDir,dirHorznIrradNormHourlyMod(:,14));
cropStructLettuce = struct();
cropStructLettuce.type = 'lettuce';   % nothing else needed here
% ---- Call uncertainty runner ----
run_uncertainty_lettuce( ...
    10, ...
    tspan, ...
    numGroundX, ...
    numGroundY, ...
    Tao_daily, ...
    rhao_daily, ...
    wind_daily, ...
    gamma_var, ...
    lambda, ...
    C_pv_front, ...
    C_pv_back, ...
    daily_front_abs_totals, ...
    daily_rear_abs_totals, ...
    C_ground_top, ...
    k_ground, ...
    ground_albedo_int, ...
    C_crop, ...
    k_crop, ...
    C_ground_inner, ...
    h_conv_inner, ...
    wavelengthDir, ...
    panelWidth, ...
    panelHeight, ...
    Pr, ...
    numPanels, ...
    sig, ...
    epsilon_pv, ...
    ri, ...
    sp, ...
    daily_ePAR_totals, ...
    target_head_weight, ...
    planting_density, ...
    dry_matter_fraction, ...
    usable_fraction, ...
    max_days, ...
    theta_m, ...
    Ca, Csl, Oa, ra, rb, Pre, ...
    CT, T0, ...
    Vcmax0, g1, g0, rjv, theta, alpha, ...
    Ws_init, Wg_init, GDD_init, ...
    c_shoot, c_root, f_root, cp, ...
    c_gr_max, Tbase, c_lar, ...
    T_sky_K_day, ...
    hours_per_day, ...
    hourNumbersSorted_VF, ...
    viewFactor_GroundToPVFront, ...
    viewFactor_GroundToPVRear, ...
    cf_matrix, ...
    cropStructLettuce, ...
    daily_totals, ...
    daily_ePAR_mole_totals,k_ext_new,...
    SLAmin_new,SLAmax_new,...
        SLAexp_new,betaT_new,betaC_new);

save('GI_Yuma_APVSS.mat','groundIrradiance_Yuma_APVSS','hourNumbersSorted_Yuma_APVSS_grIr','-v7.3');

save('PI_Yuma_APVSS.mat','panelirradiance_Yuma_APVSS','hourNumbersSorted_Yuma_APVSS_pair','-v7.3');

save('PP_Yuma_APVSS.mat','PV_Yuma_APVSS', 'hourNumbersSorted_Yuma_APVSS_PV','-v7.3');

