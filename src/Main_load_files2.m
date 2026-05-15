clear; close all;

%% ==== EDIT THESE PATHS BEFORE RUNNING =============================
% PLOTS_BASE       Full path to plots_base.mat.
% HOTMIRROR_XLSX   Hot-mirror spreadsheet for this case (case 2).
% SCRATCH_DIR      Reused (geometry-only) outputs: viewfactor_groundPV,
%                  ground_shadows, skyViewFactor, sunVectors,
%                  self_shading_direct, sky_diffuse_shading.
% SCRATCH_DIR1     Case-specific outputs: panelIrradiance, groundIrradiance,
%                  spectral_components, PVPerformance, Photosynthesis.
PLOTS_BASE     = 'plots_base.mat';
HOTMIRROR_XLSX = 'Hotmirror_T_R_case2.xlsx';
SCRATCH_DIR    = fullfile(pwd, 'scratch', 'uniform_plus');
SCRATCH_DIR1   = fullfile(pwd, 'scratch', 'uniform_minus');
%% ==================================================================

load(PLOTS_BASE)
filename_HM = HOTMIRROR_XLSX;

% Read data starting from row 2 to ignore the first row (assumed to be headers or junk)

% R_unpol_smoothed (starts from A2)
R_HM_allangles = readmatrix(filename_HM, 'Sheet', 'Reflectance', 'Range', 'A2');
T_HM_allangles = readmatrix(filename_HM, 'Sheet', 'Transmittance', 'Range', 'A2');
scratchDir = SCRATCH_DIR;
scratchDir1 = SCRATCH_DIR1;
% Define the folder name
folderName = fullfile(scratchDir, 'viewfactor_groundPV');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

folderName2 = fullfile(scratchDir, 'ground_shadows');
if ~exist(folderName2, 'dir')
    mkdir(folderName2);
end

folderName3 = fullfile(scratchDir, 'skyViewFactor');
if ~exist(folderName3, 'dir')
    mkdir(folderName3);
end

folderName4 = fullfile(scratchDir1, 'panelIrradiance');
if ~exist(folderName4, 'dir')
    mkdir(folderName4);
end

folderName5 = fullfile(scratchDir1, 'groundIrradiance');
if ~exist(folderName5, 'dir')
    mkdir(folderName5);
end

folderName6 = fullfile(scratchDir, 'sunVectors');
if ~exist(folderName6, 'dir')
    mkdir(folderName6);
end

folderName7 = fullfile(scratchDir, 'self_shading_direct');
if ~exist(folderName7, 'dir')
    mkdir(folderName7);
end

folderName8 = fullfile(scratchDir, 'sky_diffuse_shading');
if ~exist(folderName8, 'dir')
    mkdir(folderName8);
end

folderName9 = fullfile(scratchDir1, 'spectral_components');
if ~exist(folderName9, 'dir')
    mkdir(folderName9);
end

folderName10 = fullfile(scratchDir1, 'PVPerformance');
if ~exist(folderName10, 'dir')
    mkdir(folderName10);
end
folderName11 = fullfile(scratchDir, 'Photosynthesis');
if ~exist(folderName11, 'dir')
    mkdir(folderName11);
end
% cpus = str2double(getenv('SLURM_CPUS_PER_TASK'));
% if isnan(cpus) || cpus < 1
%     cpus = feature('numcores');
% end
% 
% c = parcluster('local');
% if c.NumWorkers < cpus
%     c.NumWorkers = cpus;
% end
% 
% pool = gcp('nocreate');
% if isempty(pool)
%     parpool(c, cpus);
% elseif pool.NumWorkers ~= cpus
%     delete(pool);
%     parpool(c, cpus);
% end

%% Hour range
hours = [1:744  5832:8760]';
% Keep only daytime hours
dayMask = IH(hours,1) > 0;
hours   = hours(dayMask);

nHours  = numel(hours);

%% Preallocate struct array
fileList(nHours) = struct( ...
    'sunVec', '', ...
    'groundShadow', '', ...
    'viewFactor', '', ...
    'selfShade', '', ...
    'skyDiffuse', '', ...
    'skyView', '' );

%% Build filenames ONLY for daytime hours
for k = 1:nHours
    i = hours(k);

    fileList(k).sunVec = fullfile(folderName6, ...
        sprintf('sunVector_hour_%d.mat', i));

    fileList(k).groundShadow = fullfile(folderName2, ...
        sprintf('groundshadow_hour_%d.mat', i));

    fileList(k).viewFactor = fullfile(folderName, ...
        sprintf('viewFactor_hour_%d.mat', i));

    fileList(k).selfShade = fullfile(folderName7, ...
        sprintf('self_shaded_direct_hour_%d.mat', i));

    fileList(k).skyDiffuse = fullfile(folderName8, ...
        sprintf('sky_diffuse_shading_hour_%d.mat', i));

    fileList(k).skyView = fullfile(folderName3, ...
        sprintf('skyViewFactor_hour_%d.mat', i));
end

%% Verify files exist (daytime only)

filesPerHour = 6;
missingFiles = cell(nHours * filesPerHour, 1);
idx = 0;

for k = 1:nHours
    f = struct2cell(fileList(k));

    for j = 1:numel(f)
        if ~isfile(f{j})
            idx = idx + 1;
            missingFiles{idx} = f{j};
        end
    end
end

missingFiles = missingFiles(1:idx);

if ~isempty(missingFiles)
    fprintf('\nMissing %d files (daytime only):\n\n', idx);
    fprintf('%s\n', missingFiles{:});
    error('Fix missing daytime files before running parfor.');
end

%% Simulating available light across the 
tic
dirMat = dirHorznIrradNormHourlyMod(:, hours);
difMat = difHorznIrradNormHourlyMod(:, hours);
parfor k = 1:nHours

    i = hours(k);
    
    % --- your original lines ---
    doy    = ceil(i/24);
    HExtra = pvl_extraradiation(doy);

    % Skip night hours
    if IH(i,1) <= 0
        continue
    end

    %% ---- Load Files Using Prebuilt List ----
    filename6  = fileList(k).sunVec;
    filename2  = fileList(k).groundShadow;
    filenameVF = fileList(k).viewFactor;
    filename7  = fileList(k).selfShade;
    filename8  = fileList(k).skyDiffuse;
    filename3  = fileList(k).skyView;

    % Load sun vectors
    S_sun = loadSunVectors(filename6, trackingType);

    % Load ground shadow
    S_shadow = loadGroundShadow(filename2, trackingType);

    % Load view factors
    S_vf = loadViewFactorGroundPV(filenameVF);

    % Load self shading
    S_self = loadSelfShading_Direct(filename7, trackingType);

    % Load sky diffuse reduction
    S_skyDiffuse = loadSkyDiffuseReduction(filename8);

    % Load sky view factors
    S_skyView = loadSkyViewFactor(filename3, trackingType);

    %% ---- Now use S_* variables exactly as in your original code ----
    n_panel      = S_sun.n_panel;
    sunDir       = S_sun.sunDir;
    v_shadow     = S_sun.v_shadow;
    panelTilt    = S_sun.panelTilt;
    panelAzimuth = S_sun.panelAzimuth;

    n_hm1        = S_sun.n_hm1;
    n_hm2        = S_sun.n_hm2;
    hm1Tilt      = S_sun.hm1Tilt;
    hm2Tilt      = S_sun.hm2Tilt;
    hm1Azimuth   = S_sun.hm1Azimuth;
    hm2Azimuth   = S_sun.hm2Azimuth;

    panelCenters   = S_shadow.panelCenters;
    panelCorners   = S_shadow.panelCorners;
    shadowPolygons = S_shadow.shadowPolygons;
    hm1Center      = S_shadow.hm1Center;
    hm2Center      = S_shadow.hm2Center;
    hmCorners      = S_shadow.hmCorners;
    hmShadows      = S_shadow.hmShadows;
    entityShadows  = S_shadow.entityShadows;

    viewFactor_FrontPVGround  = cellfun(@full, S_vf.viewFactor_FrontPVGround_sparse,  'UniformOutput', false);
    viewFactor_RearPVGround   = cellfun(@full, S_vf.viewFactor_RearPVGround_sparse,   'UniformOutput', false);
    viewFactor_FrontHM1Ground = cellfun(@full, S_vf.viewFactor_FrontHM1Ground_sparse, 'UniformOutput', false);
    viewFactor_RearHM1Ground  = cellfun(@full, S_vf.viewFactor_RearHM1Ground_sparse,  'UniformOutput', false);
    viewFactor_FrontHM2Ground = cellfun(@full, S_vf.viewFactor_FrontHM2Ground_sparse, 'UniformOutput', false);
    viewFactor_RearHM2Ground  = cellfun(@full, S_vf.viewFactor_RearHM2Ground_sparse,  'UniformOutput', false);

    shadingFractions_direct    = S_self.shadingFractions_direct;
    shadingFractionsHM1_direct = S_self.shadingFractionsHM1_direct;
    shadingFractionsHM2_direct = S_self.shadingFractionsHM2_direct;

    reduction_dualHM1 = S_skyDiffuse.reduction_dualHM1;
    reduction_dualHM2 = S_skyDiffuse.reduction_dualHM2;
    reduction_dual    = S_skyDiffuse.reduction_dual;
    SVF_ground        = S_skyDiffuse.SVF_ground;

    SVF_front    = S_skyView.SVF_front;
    SVF_back     = S_skyView.SVF_back;
    SVF_HM1front = S_skyView.SVF_HM1front;
    SVF_HM1back  = S_skyView.SVF_HM1back;
    SVF_HM2front = S_skyView.SVF_HM2front;
    SVF_HM2back  = S_skyView.SVF_HM2back;
        
      
        
        
        
        [AOI_HM1, AOI_HM1_Hor, R_HM1, T_HM1,A_HM1] = computeHMOpticsDualFace( ...
        hm1Tilt, hm1Azimuth, zenith(i,1), az_s_matrix(i,1), ...
        R_HM_allangles, T_HM_allangles,A_HM_allangles, wavelength_optics, wavelengthDir);

        [AOI_HM2, AOI_HM2_Hor, R_HM2, T_HM2,A_HM2] = computeHMOpticsDualFace( ...
        hm2Tilt, hm2Azimuth, zenith(i,1), az_s_matrix(i,1), ...
        R_HM_allangles, T_HM_allangles,A_HM_allangles, wavelength_optics, wavelengthDir);
    
        [AOI_Panel, AOI_Panel_Hor, R_Panel, T_Panel,A_Panel] =...
         computeHMOpticsDualFace_Panel( ...
        sunDir,n_hm1,n_hm2,n_panel,panelTilt, panelAzimuth,...
        R_stack_allangles, T_stack_allangles, ...
        A_cSi_allangles, ...
        wavelength_optics, wavelengthDir);

        AM = pvl_relativeairmass(zenith(i,1));
        AM(isnan(AM)) = 20;
        
        

        VF_HM_to_Panel = ...
        computeRayBasedHMToPanelVF(hmCorners, ...
        panelCorners,n_panel,n_hm1,n_hm2, ...
        N_rays, maxViewAngle_deg);
        
        VF_Panel_to_HM = computeReciprocalVF_PanelToHM( ...
        VF_HM_to_Panel, A_HM,A_PV);

        [viewFactorCell_GroundToHM1Front,...
            viewFactorCell_GroundToHM1Front_sparse] ...
            = convertPVToGroundToGroundToPV...
            (viewFactor_FrontHM1Ground, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, ...
        groundYmin, groundYmax, numGroundX, numGroundY);

        [viewFactorCell_GroundToHM1Rear,...
            viewFactorCell_GroundToHM1Rear_sparse] ...
            = convertPVToGroundToGroundToPV(viewFactor_RearHM1Ground, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, groundYmin, ...
        groundYmax, numGroundX, numGroundY);  

        [viewFactorCell_GroundToHM2Front,...
            viewFactorCell_GroundToHM2Front_sparse] ...
            = convertPVToGroundToGroundToPV...
            (viewFactor_FrontHM2Ground, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, ...
        groundYmin, groundYmax, numGroundX, numGroundY);

        [viewFactorCell_GroundToHM2Rear,...
            viewFactorCell_GroundToHM2Rear_sparse] ...
            = convertPVToGroundToGroundToPV(viewFactor_RearHM2Ground, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, groundYmin, ...
        groundYmax, numGroundX, numGroundY);  

        [viewFactorCell_GroundToPVFront,...
            viewFactorCell_GroundToPVFront_sparse] ...
            = convertPVToGroundToGroundToPV...
            (viewFactor_FrontPVGround, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, ...
        groundYmin, groundYmax, numGroundX, numGroundY);

        [viewFactorCell_GroundToPVRear,...
            viewFactorCell_GroundToPVRear_sparse] ...
            = convertPVToGroundToGroundToPV(viewFactor_RearPVGround, ...
        panelWidth, panelHeight, ...
        groundXmin, groundXmax, groundYmin, ...
        groundYmax, numGroundX, numGroundY);
                           

        % Construct the full file path within the subfolder
        
        % Save the variable in the specified subfolder
       
        hm1Corners = cell(size(hmCorners));  % Preallocate cell array
        hm2Corners = cell(size(hmCorners));  % Preallocate cell array
        for t = 1:length(hmCorners)
            hm1Corners{t} = hmCorners{t}{1};  % Extract the first element from the 1x2 cell
            hm2Corners{t} = hmCorners{t}{2};
        end
        
        
        
        

        

        % Reshape shadingFractions_direct into a 12x1 column vector
        shadingFractionsHM1_direct_reshaped = ...
            reshape(shadingFractionsHM1_direct', ...
            numPanels, 1);  % 12x1
         shadingFractionsHM2_direct_reshaped = ...
            reshape(shadingFractionsHM2_direct', ...
            numPanels, 1);  % 12x1
         shadingFractions_direct_reshaped = ...
            reshape(shadingFractions_direct', ...
            numPanels, 1);  % 12x1
           
       
         
        hmCorners_total = cell(numPanels, 1);

        for iteration = 1:numPanels
            hmCorners_total{iteration} = [hm1Corners{iteration}; hm2Corners{iteration}];  % Combine 4x3 + 4x3 = 8x3
        end         
        
        
        
        % unique_LAI = unique(LAI);
        % [sun_tables, shade_tables] = ...
        %     run_prosail_batch(unique_LAI(1,1), );
        
        [sun_ref_tables, shade_ref_tables] = ...
            run_prosail_batch(ref_LAIs', zenith(i,1), az_s_matrix(i,1));

        [sun_interp,shade_interp,...
        sun_albedo, shade_albedo, ...
        dirHorznIrradNormHourlyMod_new, difHorznIrradNormHourlyMod_new] = ...
        calculate_spectral_components(dirMat(:,k), ...
        difMat(:,k), ...
        wavelengthDir, wavelength,sun_ref_tables,shade_ref_tables);

        [refl_sun, trans_sun] = interpolate_canopy_spectra(LAIassump, ref_LAIs, sun_interp);
        %sig_sun = cellfun(@(a, b) a + b, refl_sun, trans_sun, 'UniformOutput', false);
        %sig_sun = refl_sun + trans_sun;

        [refl_shade, trans_shade] = interpolate_canopy_spectra(LAIassump, ref_LAIs, shade_interp);
        %sig_shade = cellfun(@(a, b) a + b, refl_shade, trans_shade, 'UniformOutput', false);
        %sig_shade = refl_shade + trans_shade;
        
        
        
        [Front_Irradiance,Front_Irradiance_HM1,Front_Irradiance_HM2,...
         Rear_Irradiance,Rear_Irradiance_HM1,Rear_Irradiance_HM2,...
         albedo_spec_per_PV,albedo_spec_per_HM1,albedo_spec_per_HM2,...
         albedo_spec_per_PV_rear,albedo_spec_per_HM1_rear,...
         albedo_spec_per_HM2_rear,out] = ...
         computePanelIrradianceBifacialRayVF_APVSS(DHI(i,1), DNI(i,1), HExtra, ...
         zenith(i,1), az_s_matrix(i,1), AM, ...
         hm1Tilt,hm2Tilt,hm1Azimuth,hm2Azimuth, ...
         AOI_HM1, ...
         R_HM1,R_HM2,A_Panel, ...
         reduction_dualHM1,reduction_dualHM2,...
         numPanelRows,...
         numPanelCols,dirHorznIrradNormHourlyMod_new, ....
         difHorznIrradNormHourlyMod_new,wavelengthDir,viewFactorCell_GroundToPVFront,...
         viewFactorCell_GroundToPVRear,refl_sun,...
         refl_shade,entityShadows,groundX,groundY,....
         shadingFractionsHM1_direct,VF_HM_to_Panel,...
         SVF_HM1front,T_HM1,R_Panel,...
         A_HM1,A_HM2,...
         viewFactorCell_GroundToHM1Front,...
         viewFactorCell_GroundToHM2Front,AOI_HM2,shadingFractionsHM2_direct,...
         viewFactorCell_GroundToHM1Rear,viewFactorCell_GroundToHM2Rear,SVF_ground);
        % Front component              
        IB_Front             = out.IB_Front;
        IB_Front_HM1         = out.IB_Front_HM1;
        IB_Front_HM2         = out.IB_Front_HM2;
        I_Alb_Front          = out.I_Alb_Front;
        I_Alb_HM1_Front      = out.I_Alb_HM1_Front;
        I_Alb_HM2_Front      = out.I_Alb_HM2_Front;
        ID_Iso_Front         = out.ID_Iso_Front;
        ID_Iso_Front_HM1     = out.ID_Iso_Front_HM1;
        ID_Iso_Front_HM2     = out.ID_Iso_Front_HM2;
        ID_Cir_Front         = out.ID_Cir_Front;
        ID_Cir_Front_HM1     = out.ID_Cir_Front_HM1;
        ID_Cir_Front_HM2     = out.ID_Cir_Front_HM2;
        ID_Hor_Front         = out.ID_Hor_Front;
        ID_Hor_Front_HM1     = out.ID_Hor_Front_HM1;
        ID_Hor_Front_HM2     = out.ID_Hor_Front_HM2;
        IB_Front_abs         = out.IB_Front_abs;
        IB_Front_abs_HM1     = out.IB_Front_abs_HM1;
        IB_Front_abs_HM2     = out.IB_Front_abs_HM2;
        ID_Iso_Front_abs     = out.ID_Iso_Front_abs;
        ID_Iso_Front_abs_HM1 = out.ID_Iso_Front_abs_HM1;
        ID_Iso_Front_abs_HM2 = out.ID_Iso_Front_abs_HM2;
        ID_Cir_Front_abs     = out.ID_Cir_Front_abs;
        ID_Cir_Front_abs_HM1 = out.ID_Cir_Front_abs_HM1;
        ID_Cir_Front_abs_HM2 = out.ID_Cir_Front_abs_HM2;
        ID_Hor_Front_abs     = out.ID_Hor_Front_abs;
        ID_Hor_Front_abs_HM1 = out.ID_Hor_Front_abs_HM1;
        ID_Hor_Front_abs_HM2 = out.ID_Hor_Front_abs_HM2;
        albedo_spec_per_PV_abs      = out.albedo_spec_per_PV_abs;
        albedo_spec_per_HM1_abs     = out.albedo_spec_per_HM1_abs;
        albedo_spec_per_HM2_abs     = out.albedo_spec_per_HM2_abs;
        albedo_per_PV_abs           = out.albedo_per_PV_abs;
        albedo_per_HM1_Front_abs          = out.albedo_per_HM1_Front_abs;
        albedo_per_HM2_Front_abs          = out.albedo_per_HM2_Front_abs;
        I_Alb_Front_abs             = out.I_Alb_Front_abs;
        I_Alb_HM1_Front_abs         = out.I_Alb_HM1_Front_abs;
        I_Alb_HM2_Front_abs         = out.I_Alb_HM2_Front_abs;
        Front_Irradiance_abs        = out.Front_Irradiance_abs;
        Front_Irradiance_HM1_abs    = out.Front_Irradiance_HM1_abs;
        Front_Irradiance_HM2_abs    = out.Front_Irradiance_HM2_abs;

        %Rear component
        IB_Back              = out.IB_Back;
        IB_Back_HM1          = out.IB_Back_HM1;
        IB_Back_HM2          = out.IB_Back_HM2;
        I_Alb_Rear           = out.I_Alb_Rear;
        I_Alb_HM1_Rear       = out.I_Alb_HM1_Rear;
        I_Alb_HM2_Rear       = out.I_Alb_HM2_Rear;
        ID_Iso_Rear          = out.ID_Iso_Rear;
        ID_Iso_Rear_HM1      = out.ID_Iso_Rear_HM1;
        ID_Iso_Rear_HM2      = out.ID_Iso_Rear_HM2;
        ID_Cir_Rear          = out.ID_Cir_Rear;
        ID_Cir_Rear_HM1      = out.ID_Cir_Rear_HM1;
        ID_Cir_Rear_HM2      = out.ID_Cir_Rear_HM2;
        ID_Hor_Rear          = out.ID_Hor_Rear;
        ID_Hor_Rear_HM1      = out.ID_Hor_Rear_HM1;
        ID_Hor_Rear_HM2      = out.ID_Hor_Rear_HM2;
        IB_Back_abs          = out.IB_Back_abs;
        IB_Back_HM1_abs      = out.IB_Back_HM1_abs;
        IB_Back_HM2_abs      = out.IB_Back_HM2_abs;
        ID_Iso_Rear_abs      = out.ID_Iso_Rear_abs;
        ID_Iso_Rear_abs_HM1  = out.ID_Iso_Rear_abs_HM1;
        ID_Iso_Rear_abs_HM2  = out.ID_Iso_Rear_abs_HM2;
        ID_Cir_Rear_abs      = out.ID_Cir_Rear_abs;
        ID_Cir_Rear_abs_HM1  = out.ID_Cir_Rear_abs_HM1;
        ID_Cir_Rear_abs_HM2  = out.ID_Cir_Rear_abs_HM2;
        ID_Hor_Rear_abs      = out.ID_Hor_Rear_abs;
        ID_Hor_Rear_abs_HM1  = out.ID_Hor_Rear_abs_HM1;
        ID_Hor_Rear_abs_HM2  = out.ID_Hor_Rear_abs_HM2;
        I_Alb_Rear_abs             = out.I_Alb_Rear_abs;
        I_Alb_HM1_Rear_abs         = out.I_Alb_HM1_Rear_abs;
        I_Alb_HM2_Rear_abs         = out.I_Alb_HM2_Rear_abs;
        albedo_spec_per_PV_rear_abs   = out.albedo_spec_per_PV_rear_abs;
        albedo_spec_per_HM1_rear_abs  = out.albedo_spec_per_HM1_rear_abs;
        albedo_spec_per_HM2_rear_abs  = out.albedo_spec_per_HM2_rear_abs;
        albedo_per_PV_rear_abs        = out.albedo_per_PV_rear_abs;
        albedo_per_HM1_rear_abs       = out.albedo_per_HM1_rear_abs;
        albedo_per_HM2_rear_abs       = out.albedo_per_HM2_rear_abs;
        Rear_Irradiance_abs           = out.Back_Irradiance_abs;
        Rear_Irradiance_HM1_abs       = out.Rear_Irradiance_HM1_abs;
        Rear_Irradiance_HM2_abs       = out.Rear_Irradiance_HM2_abs;

        
        [IB_Front_inten,ID_Iso_Front_inten,ID_Cir_Front_inten,... 
        ID_Hor_Front_inten,IB_Back_inten,ID_Iso_Rear_inten,...
        ID_Cir_Rear_inten,ID_Hor_Rear_inten,...
        IB_HM1_Front_inten,ID_HM1_Iso_Front_inten,ID_HM1_Cir_Front_inten,... 
        ID_HM1_Hor_Front_inten,IB_HM1_Back_inten,ID_HM1_Iso_Rear_inten,...
        ID_HM1_Cir_Rear_inten,ID_HM1_Hor_Rear_inten,...
        IB_HM2_Front_inten,ID_HM2_Iso_Front_inten,ID_HM2_Cir_Front_inten,... 
        ID_HM2_Hor_Front_inten,IB_HM2_Back_inten,ID_HM2_Iso_Rear_inten,...
        ID_HM2_Cir_Rear_inten,ID_HM2_Hor_Rear_inten] = ...
        computePanelHMIntensity(numPanels,IB_Front,IB_Front_HM1,IB_Front_HM2,...
        ID_Iso_Front,...
        ID_Iso_Front_HM1,ID_Iso_Front_HM2,...
        ID_Cir_Front,ID_Cir_Front_HM1,...
        ID_Cir_Front_HM2,ID_Hor_Front,...
        ID_Hor_Front_HM1,ID_Hor_Front_HM2,...
        IB_Back,IB_Back_HM1,...
        IB_Back_HM2,ID_Iso_Rear,...
        ID_Iso_Rear_HM1,ID_Iso_Rear_HM2,...
        ID_Cir_Rear,ID_Cir_Rear_HM1,...
        ID_Cir_Rear_HM2,...
        ID_Hor_Rear,ID_Hor_Rear_HM1,...
        ID_Hor_Rear_HM2,wavelengthDir);

        
        filename9 = fullfile(folderName9, sprintf('spectralComponents_hour_%d.mat', i));
        
        saveSpectralFactors_APVSS(filename9,R_HM1, R_HM2, R_Panel,A_HM1,A_HM2,...
        A_Panel,T_HM1,T_HM2,T_Panel,...
        dirHorznIrradNormHourlyMod_new, difHorznIrradNormHourlyMod_new);
        

        %Construct the full file path within the subfolder
        [panel_reflected_front, panel_reflected_rear] = ...
            computePanelReflectedSpectral( ...
        numPanels, groundX,Front_Irradiance,...
        Rear_Irradiance,viewFactor_FrontPVGround,...
        viewFactor_RearPVGround,wavelengthDir,...
        R_Panel.Front.Albedo,...
        R_Panel.Back.Albedo);

        [hm1_reflected_front, hm1_reflected_rear] = ...
            computePanelReflectedSpectral( ...
        numPanels, groundX,Front_Irradiance_HM1,...
        Rear_Irradiance_HM1,viewFactor_FrontHM1Ground,...
        viewFactor_RearHM1Ground,wavelengthDir,...
        R_HM1.Front.Albedo,...
        R_HM1.Back.Albedo);

        [hm2_reflected_front, hm2_reflected_rear] = ...
            computePanelReflectedSpectral( ...
        numPanels, groundX,Front_Irradiance_HM2,...
        Rear_Irradiance_HM2,viewFactor_FrontHM2Ground,...
        viewFactor_RearHM2Ground,wavelengthDir,...
        R_HM2.Front.Albedo,...
        R_HM2.Back.Albedo);
        
        [Front_Irradiance_inten, Rear_Irradiance_inten, ...
            panel_reflected_front_inten, panel_reflected_rear_inten] = ...
        computeTotalPanelIntensity(Front_Irradiance, Rear_Irradiance, ...
        panel_reflected_front, panel_reflected_rear, wavelengthDir);

        [Front_Irradiance_inten_HM1, Rear_Irradiance_inten_HM1, ...
            HM1_reflected_front_inten, HM1_reflected_rear_inten] = ...
        computeTotalPanelIntensity(Front_Irradiance_HM1, Rear_Irradiance_HM1, ...
        hm1_reflected_front, hm1_reflected_rear, wavelengthDir);

        [Front_Irradiance_inten_HM2, Rear_Irradiance_inten_HM2, ...
            HM2_reflected_front_inten, HM2_reflected_rear_inten] = ...
        computeTotalPanelIntensity(Front_Irradiance_HM2, Rear_Irradiance_HM2, ...
        hm2_reflected_front, hm2_reflected_rear, wavelengthDir);

        filename4 = fullfile(folderName4, sprintf('panelIrradiance_hour_%d.mat', i));

       % Save the variable in the specified subfolder
        savePanelHMIrradiance(filename4,Front_Irradiance_inten, ...
            Front_Irradiance_inten_HM1,Rear_Irradiance_inten_HM1, Rear_Irradiance_inten,...
            Front_Irradiance_inten_HM2,...
            Rear_Irradiance_inten_HM2,panel_reflected_front_inten, ...
            HM1_reflected_front_inten, ...
            HM1_reflected_rear_inten,...
            panel_reflected_rear_inten,...
            HM2_reflected_front_inten, ...
            HM2_reflected_rear_inten,...
            IB_Front_inten,...
            IB_HM1_Front_inten,IB_HM1_Back_inten,...
            I_Alb_Front,I_Alb_HM1_Front,...
            I_Alb_HM1_Rear,ID_Iso_Front_inten,...
            ID_HM1_Iso_Front_inten, ID_HM1_Iso_Rear_inten,...
            ID_Cir_Front_inten,...
            ID_HM1_Cir_Front_inten, ID_HM1_Cir_Rear_inten,...
            ID_Hor_Front_inten,...
            ID_HM1_Hor_Front_inten, ID_HM1_Hor_Rear_inten,...
            IB_Back_inten,...
            IB_HM2_Front_inten,IB_HM2_Back_inten,...
            I_Alb_Rear,I_Alb_HM2_Front,...
            I_Alb_HM2_Rear,ID_Iso_Rear_inten,...
            ID_HM2_Iso_Front_inten, ID_HM2_Iso_Rear_inten,...
            ID_Cir_Rear_inten, ...
            ID_HM2_Cir_Front_inten, ID_HM2_Cir_Rear_inten, ...
            ID_Hor_Rear_inten,...
            ID_HM2_Hor_Front_inten, ID_HM2_Hor_Rear_inten);
              
        panelResults = ...
        computePVPanelPerformance(...
        wavelengthDir, E0, ModuleParameters, ...
        wind_vel(i,1), Tao(i,1), numPanels, EgRef, dEgdt,Front_Irradiance_abs,bifaciality_factor,Rear_Irradiance_abs,bifacialModuleNamesSorted);

        filename10 = fullfile(folderName10, sprintf('PVPanelPerformance_hour_%d.mat', i));
       % Save the variable in the specified subfolder
        savePVPanelPerformance(filename10,panelResults,Front_Irradiance_abs,Rear_Irradiance_abs);
        % Convert angles to radians for calculations
        zenithRad = deg2rad(zenith(i,1));
        
        [~, ID_Iso_sky, ID_Cir_sky, ID_Hor_sky] = ...
        pvl_perez(0, 0, ...
        DHI(i,1), DNI(i,1), HExtra, zenith(i,1), az_s_matrix(i,1), AM,...
        difHorznIrradNormHourlyMod_new,dirHorznIrradNormHourlyMod_new,...
        wavelengthDir, '1990');
        

       [PAR_inten, ePAR_inten, IR_inten,Ri_inten, PAR_inten_mmole, ...
           ePAR_inten_mmole,Ri_spec,shadingMaskhm1,...
        shadingMaskhm2,directSolar,diffuseSolar,idx_PAR,idx_ePAR] = ...
        computeGroundIrradiance_APVSS(...
        groundX, groundY, DNI(i,1),IB_Front_HM1, IB_Front_HM2, ID_Iso_Front_HM1, ID_Iso_Front_HM2,...
        ID_Cir_Front_HM1, ID_Cir_Front_HM2,...
        ID_Hor_Front_HM1, ID_Hor_Front_HM2, ...
        zenithRad, hmShadows, ...
        viewFactorCell_GroundToHM1Rear, viewFactorCell_GroundToHM2Rear, ...
        panel_reflected_front, panel_reflected_rear, ...
        dirHorznIrradNormHourlyMod_new, difHorznIrradNormHourlyMod_new, wavelengthDir,...
        T_HM1,T_HM2,ID_Iso_sky, ID_Cir_sky, ID_Hor_sky,SVF_ground);
        
        filename5 = fullfile(folderName5, sprintf('irradianceGround_hour_%d.mat', i)); 
        saveGroundIrradiance_APVSS(filename5,PAR_inten, ePAR_inten,...
            IR_inten,Ri_inten,PAR_inten_mmole, ePAR_inten_mmole,...
            shadingMaskhm1,...
            shadingMaskhm2);

        % resultsPhotosynthesis = runPhotosynthesisGrid(Ri_spec, Tao(i,1), ...
        % ri, qy, ab, Ca, Csl, ra, rb, Pre, RH, CT, T0, Vcmax0, ...
        % Oa, g1, g0, rjv, theta, alpha,wind_vel(i,1),...
        % LAIassump);
        % 
        % filename11 = fullfile(folderName11, sprintf('photosynthesis_hour_%d.mat', i));
        % savePhotosynthesis(filename11, resultsPhotosynthesis);

end


toc