function [Front_Irradiance,Front_Irradiance_HM1,Front_Irradiance_HM2,...
    Rear_Irradiance,Rear_Irradiance_HM1,Rear_Irradiance_HM2,...
    albedo_spec_per_PV,albedo_spec_per_HM1,albedo_spec_per_HM2,...
    albedo_spec_per_PV_rear,albedo_spec_per_HM1_rear,...
    albedo_spec_per_HM2_rear,out] = ...
    computePanelIrradianceBifacialRayVF_APVSS(DHI, DNI, HExtra, zenith, azimuth, AM, ...
    hm1Tilt,hm2Tilt,hm1Azimuth,hm2Azimuth, ...
    AOI_HM1, ...
    R_HM1,R_HM2,A_Panel, ...
    reduction_dualHM1,reduction_dualHM2,...
    numPanelRows,...
    numPanelColumns,dirHorznIrradNormHourlyMod, ....
    difHorznIrradNormHourlyMod,wavelengthDir,...
    viewFactorCell_GroundToPVFront,...
    viewFactorCell_GroundToPVRear,refl_sun,...
        refl_shade,entityShadows,groundX,groundY,....
        shadingFractionsHM1_direct,VF_HM_to_Panel,SVF_HM1front,...
        T_HM1,R_Panel,A_HM1,A_HM2,...
        viewFactorCell_GroundToHM1Front,...
        viewFactorCell_GroundToHM2Front,AOI_HM2,shadingFractionsHM2_direct,...
        viewFactorCell_GroundToHM1Rear,viewFactorCell_GroundToHM2Rear,SVF_ground)



        % Determine reduction factor for diffuse shading
        %reduction_dif = reduction_dual;
    
        % Preallocate cell array
        IB_Front = cell(numPanelRows*numPanelColumns, 1);
        IB_Front_abs = cell(numPanelRows*numPanelColumns, 1);
        IB_Front_HM1 = cell(numPanelRows*numPanelColumns, 1);
        IB_Front_abs_HM1 = cell(numPanelRows*numPanelColumns, 1);
        
        IB_Back_HM1 = cell(numPanelRows*numPanelColumns, 1);
        IB_Back_HM1_abs = cell(numPanelRows*numPanelColumns, 1);

        for p = 1:numPanelRows*numPanelColumns
            
            IB_Front{p} = DNI.* cosd(AOI_HM1.Front) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  R_HM1.Front.Beam.*dirHorznIrradNormHourlyMod;
            IB_Front_abs{p} = DNI.* cosd(AOI_HM1.Front) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  A_Panel.Front.Beam.* ...
                  R_HM1.Front.Beam.*dirHorznIrradNormHourlyMod;
            IB_Front_HM1{p} = DNI.* cosd(AOI_HM1.Front) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  dirHorznIrradNormHourlyMod;
            IB_Front_abs_HM1{p} = DNI.* cosd(AOI_HM1.Front) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  A_HM1.Front.Beam.* ...
                  dirHorznIrradNormHourlyMod;
             IB_Back_HM1{p} = DNI.* cosd(AOI_HM1.Back) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  dirHorznIrradNormHourlyMod;
            IB_Back_HM1_abs{p} = DNI.* cosd(AOI_HM1.Back) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  A_HM1.Back.Beam.* ...
                 dirHorznIrradNormHourlyMod;

            IB_Back_HM1{p}(IB_Back_HM1{p} < 0) = 0;    
            IB_Back_HM1_abs{p}(IB_Back_HM1_abs{p} < 0) = 0;
        end

% --- Sky Diffuse (Front)
[~, Iso_Front, Cir_Front, Hor_Front] = ...
    pvl_perez(hm1Tilt, hm1Azimuth, ...
    DHI, DNI, HExtra, zenith, azimuth, AM,...
    difHorznIrradNormHourlyMod,dirHorznIrradNormHourlyMod,...
    wavelengthDir, '1990');
[~, Iso_Rear, Cir_Rear, Hor_Rear] = ...
    pvl_perez(180-hm1Tilt, mod(hm1Azimuth + 180, 360), ...
    DHI, DNI, HExtra, zenith, azimuth, AM,...
    difHorznIrradNormHourlyMod,dirHorznIrradNormHourlyMod,...
    wavelengthDir, '1990');
VFHMPanel = VF_HM_to_Panel{1,1};
ID_Iso_Front = Iso_Front .* ...
               R_HM1.Front.Iso .* ...
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Cir_Front = Cir_Front .* ...
               R_HM1.Front.Beam .* ...
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Hor_Front = Hor_Front .* ...
               R_HM1.Front.Hor .* ...
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Iso_Front_abs = Iso_Front .* ...
               R_HM1.Front.Iso .* ...
               reduction_dualHM1 .*A_Panel.Front.Iso.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Cir_Front_abs = Cir_Front .* ...
               R_HM1.Front.Beam .* ...
               reduction_dualHM1 .*A_Panel.Front.Beam.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);


ID_Hor_Front_abs = Hor_Front .* ...
               R_HM1.Front.Hor .* ...
               reduction_dualHM1 .*A_Panel.Front.Hor.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);



ID_Iso_Front_HM1 = Iso_Front .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Front_HM1 = Cir_Front .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Front_HM1 = Hor_Front .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Iso_Front_abs_HM1 = Iso_Front .* ...
               reduction_dualHM1 .*A_HM1.Front.Iso.* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Front_abs_HM1 = Cir_Front .* ...
               reduction_dualHM1 .*A_HM1.Front.Beam.* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Front_abs_HM1 = Hor_Front .* ...
               reduction_dualHM1 .*A_HM1.Front.Hor.* ...
               difHorznIrradNormHourlyMod;

ID_Iso_Rear_HM1 = Iso_Rear .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Rear_HM1 = Cir_Rear .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Rear_HM1 = Hor_Rear .* ...               
               reduction_dualHM1 .* ...
               difHorznIrradNormHourlyMod;
ID_Iso_Rear_abs_HM1 = Iso_Rear .* ...
               reduction_dualHM1 .*A_HM1.Back.Iso.* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Rear_abs_HM1 = Cir_Rear .* ...
               reduction_dualHM1 .*A_HM1.Back.Beam.* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Rear_abs_HM1 = Hor_Rear .* ...
               reduction_dualHM1 .*A_HM1.Back.Hor.* ...
               difHorznIrradNormHourlyMod;
% --- Ray-Traced View Factors


% --- Albedo via ViewFactorCell + sun/shade reflectance ---

% ====== Advanced Spectral Albedo Calculation (Front) ======
nBands = length(wavelengthDir);
numPanels = length(viewFactorCell_GroundToPVFront);
[numGroundY, numGroundX] = size(viewFactorCell_GroundToPVFront{1});

% --- Step 1: Determine shaded areas on the ground ---
is_shaded = false(numGroundY, numGroundX);
for s = 1:length(entityShadows)
    poly = entityShadows{s};
    in = inpolygon(groundX, groundY, poly(:,1), poly(:,2));
    is_shaded = is_shaded | in;
end

% --- Step 2: Create reflectance maps (groundY x groundX x bands) ---
canopy_reflectance_map = zeros(numGroundY, numGroundX, nBands);

for i = 1:numGroundY
    for j = 1:numGroundX
        %LAI_val = LAI_map(i, j);
        %[~, lai_idx] = min(abs(unique_LAI - LAI_val));

        if is_shaded(i, j)
            R = refl_shade{i,j};  % nBands x 1
        else
            R = refl_sun{i,j};
        end

        canopy_reflectance_map(i, j, :) = R(:);  % broadcast into 3D array
    end
end

% --- Step 3: Compute spectrally resolved albedo for each panel ---
albedo_spec_per_PV = cell(numPanels, 1);  % each cell is nBands x 1
albedo_spec_per_PV_abs = cell(numPanels, 1);  % each cell is nBands x 1
albedo_per_PV = zeros(numPanels, 1);      % scalar per panel
albedo_per_PV_abs = zeros(numPanels, 1);      % scalar per panel

albedo_spec_per_HM1 = cell(numPanels, 1);  % each cell is nBands x 1
albedo_spec_per_HM1_abs = cell(numPanels, 1);  % each cell is nBands x 1
albedo_per_HM1 = zeros(numPanels, 1);      % scalar per panel
albedo_per_HM1_abs = zeros(numPanels, 1);      % scalar per panel

albedo_spec_per_HM1_rear = cell(numPanels, 1);  % each cell is nBands x 1
albedo_spec_per_HM1_rear_abs = cell(numPanels, 1);  % each cell is nBands x 1
albedo_per_HM1_rear = zeros(numPanels, 1);      % scalar per panel
albedo_per_HM1_rear_abs = zeros(numPanels, 1);      % scalar per panel


% Initialize irradiance map (200x200xBands)
Irradiance_map = zeros(numGroundY, numGroundX, nBands);
% Initialize the sum with zeros of the same size
% Loop through each cell and sum the matrices
for b = 1:nBands
    for i=1:numGroundY
        for j=1:numGroundX
            I_shade = ((DHI*difHorznIrradNormHourlyMod(b)).*T_HM1.Front.Iso(b,1)*...
            (1-SVF_ground(i,j)))+...
            (DHI*difHorznIrradNormHourlyMod(b)*(SVF_ground(i,j)))+...
            (DNI.* cosd(AOI_HM1.Front) .* ...
                  T_HM1.Front.Beam(b,1).*dirHorznIrradNormHourlyMod(b));
            % diffuse + direct (since hot mirrors can transmit light)
            I_sun = (DHI*difHorznIrradNormHourlyMod(b)) +...
            DNI * cosd(zenith) * dirHorznIrradNormHourlyMod(b);     % diffuse + direct
            if is_shaded(i, j)
                Irradiance_map(i,j,b) = I_shade;  % nBands x 1
            else
                Irradiance_map(i,j,b) = I_sun;
            end
    
        end
    end
end
% for p = 1:numPanels
%     VF = viewFactorCell_GroundToPV{p};  % Ny x Nx
%     % Multiply elementwise and integrate over ground
%     Irr_total = sum(sum(VF .* I_ground)) * A_cell;
%     Irr_rear_per_panel(p) = Irr_total;
% end
for k = 1:numPanels
    VF_k = viewFactorCell_GroundToPVFront{k};  % 200x200
    VF_k_HM1 = viewFactorCell_GroundToHM1Front{k};
    VF_k_HM1_rear = viewFactorCell_GroundToHM1Rear{k};
    spec_flux = zeros(nBands, 1);
    spec_flux_HM1 = zeros(nBands, 1);
    spec_flux_HM1_rear = zeros(nBands, 1);

    for b = 1:nBands
        refl_map = canopy_reflectance_map(:,:,b);       % 200x200
        irrad_b_map = Irradiance_map(:,:,b);            % 200x200
        flux_map_b = refl_map .* irrad_b_map .* VF_k;   % W contribution to panel from each pixel
        flux_map_b_HM1 = refl_map .* irrad_b_map .* VF_k_HM1;   % W contribution to panel from each pixel
        flux_map_b_HM1_rear = refl_map .* irrad_b_map .* VF_k_HM1_rear;
        spec_flux(b) = sum(flux_map_b, 'all');          % Total W for band b
        spec_flux_HM1(b) = sum(flux_map_b_HM1, 'all');
        spec_flux_HM1_rear(b) = sum(flux_map_b_HM1_rear, 'all');
    end

    % Apply albedo loss (optical loss)
    albedo_spec_per_PV{k} = ...
        spec_flux;
    albedo_spec_per_PV_abs{k} = ...
        spec_flux .*A_Panel.Front.Albedo;
    % Normalize by total view factor (acts as integration weight across the ground area)
    total_view_factor = sum(VF_k(:));

    % Apply albedo loss (optical loss)
    albedo_spec_per_HM1{k} = spec_flux_HM1;
    albedo_spec_per_HM1_abs{k} = ...
        spec_flux_HM1 .*A_HM1.Front.Albedo;
    % Normalize by total view factor (acts as integration weight across the ground area)
    total_view_factor_HM1 = sum(VF_k_HM1(:));
    
    % Apply albedo loss (optical loss)
    albedo_spec_per_HM1_rear{k} = spec_flux_HM1_rear;
    albedo_spec_per_HM1_rear_abs{k} = ...
        spec_flux_HM1_rear .*A_HM1.Back.Albedo;
    % Normalize by total view factor (acts as integration weight across the ground area)
    total_view_factor_HM1_rear = sum(VF_k_HM1_rear(:));

    if total_view_factor > 0
        albedo_per_PV(k) = trapz(wavelengthDir, albedo_spec_per_PV{k}) / total_view_factor;
        albedo_per_PV_abs(k) = trapz(wavelengthDir, albedo_spec_per_PV_abs{k}) / total_view_factor;
        albedo_spec_per_PV{k} = albedo_spec_per_PV{k} / total_view_factor;
        albedo_spec_per_PV_abs{k} = albedo_spec_per_PV_abs{k} / total_view_factor;
    else
        %warning(['Panel ', num2str(k), ': total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_PV(k) = 0;
        albedo_per_PV_abs(k) = 0;
        albedo_spec_per_PV{k} = zeros(nBands, 1);
        albedo_spec_per_PV_abs{k} = zeros(nBands, 1);
    end

    if total_view_factor_HM1 > 0
        albedo_per_HM1(k) = trapz(wavelengthDir, albedo_spec_per_HM1{k}) / total_view_factor_HM1;
        albedo_per_HM1_abs(k) = trapz(wavelengthDir, albedo_spec_per_HM1_abs{k}) / total_view_factor_HM1;
        albedo_spec_per_HM1{k} = albedo_spec_per_HM1{k} / total_view_factor_HM1;
        albedo_spec_per_HM1_abs{k} = albedo_spec_per_HM1_abs{k} / total_view_factor_HM1;
    else
        %warning(['HM1 ', num2str(k), ': total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_HM1(k) = 0;
        albedo_per_HM1_abs(k) = 0;
        albedo_spec_per_HM1{k} = zeros(nBands, 1);
        albedo_spec_per_HM1_abs{k} = zeros(nBands, 1);
    end

    if total_view_factor_HM1_rear > 0
        albedo_per_HM1_rear(k) = trapz(wavelengthDir, albedo_spec_per_HM1_rear{k}) / total_view_factor_HM1_rear;
        albedo_per_HM1_rear_abs(k) = trapz(wavelengthDir, albedo_spec_per_HM1_rear_abs{k}) / total_view_factor_HM1_rear;
        albedo_spec_per_HM1_rear{k} = albedo_spec_per_HM1_rear{k} / total_view_factor_HM1_rear;
        albedo_spec_per_HM1_rear_abs{k} = albedo_spec_per_HM1_rear_abs{k} / total_view_factor_HM1_rear;
    else
        %warning(['HM1_rear ', num2str(k), ': total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_HM1_rear(k) = 0;
        albedo_per_HM1_rear_abs(k) = 0;
        albedo_spec_per_HM1_rear{k} = zeros(nBands, 1);
        albedo_spec_per_HM1_rear_abs{k} = zeros(nBands, 1);
    end

end


% ====== Spectral Albedo Calculation (Rear) ======
albedo_spec_per_PV_rear = cell(numPanels, 1);  % nBands x 1 per panel
albedo_spec_per_PV_rear_abs = cell(numPanels, 1);  % nBands x 1 per panel
albedo_per_PV_rear = zeros(numPanels, 1);      % scalar per panel
albedo_per_PV_rear_abs = zeros(numPanels, 1);      % scalar per panel

albedo_spec_per_HM2 = cell(numPanels, 1);  % nBands x 1 per panel
albedo_spec_per_HM2_abs = cell(numPanels, 1);  % nBands x 1 per panel
albedo_per_HM2 = zeros(numPanels, 1);      % scalar per panel
albedo_per_HM2_abs = zeros(numPanels, 1);      % scalar per panel

albedo_spec_per_HM2_rear = cell(numPanels, 1);  % nBands x 1 per panel
albedo_spec_per_HM2_rear_abs = cell(numPanels, 1);  % nBands x 1 per panel
albedo_per_HM2_rear = zeros(numPanels, 1);      % scalar per panel
albedo_per_HM2_rear_abs = zeros(numPanels, 1);      % scalar per panel


for k = 1:numPanels
    VF_k_rear = viewFactorCell_GroundToPVRear{k};  % 200x200
    VF_k_HM2 = viewFactorCell_GroundToHM2Front{k};
    VF_k_HM2_rear = viewFactorCell_GroundToHM2Rear{k};
    spec_flux_rear = zeros(nBands, 1);
    spec_flux_HM2 = zeros(nBands, 1);
    spec_flux_HM2_rear = zeros(nBands, 1);

    for b = 1:nBands
        refl_map = canopy_reflectance_map(:,:,b);       % 200x200
        irrad_b_map = Irradiance_map(:,:,b);            % 200x200
        flux_map_b_rear = refl_map .* irrad_b_map .* VF_k_rear;   % W contribution to panel from each pixel
        flux_map_HM2 = refl_map .* irrad_b_map .* VF_k_HM2;   % W contribution to panel from each pixel
        flux_map_HM2_rear = refl_map .* irrad_b_map .* VF_k_HM2_rear;
        spec_flux_rear(b) = sum(flux_map_b_rear, 'all');          % Total W for band b
        spec_flux_HM2(b) = sum(flux_map_HM2, 'all');          % Total W for band b
        spec_flux_HM2_rear(b) = sum(flux_map_HM2_rear, 'all');
        
    end

    

    % Apply albedo loss (optical loss)
    albedo_spec_per_PV_rear{k} = spec_flux_rear;
    albedo_spec_per_PV_rear_abs{k} = spec_flux_rear .* A_Panel.Back.Albedo;

    albedo_spec_per_HM2{k} = spec_flux_HM2;
    albedo_spec_per_HM2_abs{k} = spec_flux_HM2 .* A_HM2.Front.Albedo;

    albedo_spec_per_HM2_rear{k} = spec_flux_HM2_rear;
    albedo_spec_per_HM2_rear_abs{k} = spec_flux_HM2_rear .* A_HM2.Back.Albedo;

    % Normalize by total view factor (acts as integration weight across the ground area)
    total_view_factor_rear = sum(VF_k_rear(:));  % unitless, scales area-wise
    total_view_factor_HM2 = sum(VF_k_HM2(:));  % unitless, scales area-wise
    total_view_factor_HM2_rear = sum(VF_k_HM2_rear(:));
    if total_view_factor_rear > 0
        albedo_per_PV_rear(k) = trapz(wavelengthDir, albedo_spec_per_PV_rear{k}) / total_view_factor_rear;
        albedo_per_PV_rear_abs(k) = trapz(wavelengthDir, albedo_spec_per_PV_rear_abs{k}) / total_view_factor_rear;
        albedo_spec_per_PV_rear{k} = albedo_spec_per_PV_rear{k} / total_view_factor_rear;
        albedo_spec_per_PV_rear_abs{k} = albedo_spec_per_PV_rear_abs{k} / total_view_factor_rear;
    else
        %warning(['Panel ', num2str(k), ' (rear): total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_PV_rear(k) = 0;
        albedo_per_PV_rear_abs(k) = 0;
        albedo_spec_per_PV_rear{k} = zeros(nBands, 1);
        albedo_spec_per_PV_rear_abs{k} = zeros(nBands, 1);
    end
    if total_view_factor_HM2 > 0
        albedo_per_HM2(k) = trapz(wavelengthDir, albedo_spec_per_HM2{k}) / total_view_factor_HM2;
        albedo_per_HM2_abs(k) = trapz(wavelengthDir, albedo_spec_per_HM2_abs{k}) / total_view_factor_HM2;
        albedo_spec_per_HM2{k} = albedo_spec_per_HM2{k} / total_view_factor_HM2;
        albedo_spec_per_HM2_abs{k} = albedo_spec_per_HM2_abs{k} / total_view_factor_HM2;
    else
        %warning(['HM2 ', num2str(k), ' (rear): total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_HM2(k) = 0;
        albedo_per_HM2_abs(k) = 0;
        albedo_spec_per_HM2{k} = zeros(nBands, 1);
        albedo_spec_per_HM2_abs{k} = zeros(nBands, 1);
    end

    if total_view_factor_HM2_rear > 0
        albedo_per_HM2_rear(k) = trapz(wavelengthDir, albedo_spec_per_HM2_rear{k}) / total_view_factor_HM2_rear;
        albedo_per_HM2_rear_abs(k) = trapz(wavelengthDir, albedo_spec_per_HM2_abs{k}) / total_view_factor_HM2_rear;
        albedo_spec_per_HM2_rear{k} = albedo_spec_per_HM2{k} / total_view_factor_HM2_rear;
        albedo_spec_per_HM2_rear_abs{k} = albedo_spec_per_HM2_abs{k} / total_view_factor_HM2_rear;
    else
        %warning(['HM2_rear ', num2str(k), ' (rear): total_view_factor = 0. Setting albedo terms to zero.']);
        albedo_per_HM2_rear(k) = 0;
        albedo_per_HM2_rear_abs(k) = 0;
        albedo_spec_per_HM2_rear{k} = zeros(nBands, 1);
        albedo_spec_per_HM2_rear_abs{k} = zeros(nBands, 1);
    end

end


% Reshape to 4×3 matrix (row × col layout)
%I_Alb_Front = reshape(I_Alb_Front, numPanelColumns, numPanelRows)';  % (4×3)

% --- Front Albedo
% I_Alb_Front = pvl_Purdue_albedo_model(panelTilt, panelAzimuth, ...
%     panelElevation / panelHeight, groundAlbedo, ...
%     DHI, DNI, HExtra, zenith, azimuth, AM, ...
%     shadingFractions_direct, reduction_dif, VF_ground, VF_alley,difHorznIrradNormHourlyMod,...
%     dirHorznIrradNormHourlyMod,wavelengthDir, '1990');
I_Alb_Front = albedo_per_PV;
I_Alb_Front(I_Alb_Front < 0) = 0;

I_Alb_Front_abs = albedo_per_PV_abs;
I_Alb_Front_abs(I_Alb_Front_abs < 0) = 0;

I_Alb_HM1_Front = albedo_per_HM1;
I_Alb_HM1_Front(I_Alb_HM1_Front < 0) = 0;
I_Alb_HM1_Front_abs = albedo_per_HM1_abs;
I_Alb_HM1_Front_abs(I_Alb_HM1_Front_abs < 0) = 0;

I_Alb_HM1_Rear = albedo_per_HM1;
I_Alb_HM1_Rear(I_Alb_HM1_Rear < 0) = 0;
I_Alb_HM1_Rear_abs = albedo_per_HM1_abs;
I_Alb_HM1_Rear_abs(I_Alb_HM1_Rear_abs < 0) = 0;



%I_Alb_Front = reshape(I_Alb_Front, numPanelColumns, numPanelRows)';  % 4×3 matrix (rows x cols)

Front_Irradiance = cell(numPanels,1);
Front_Irradiance_abs = cell(numPanels,1);

Front_Irradiance_HM1 = cell(numPanels,1);
Front_Irradiance_HM1_abs = cell(numPanels,1);

Rear_Irradiance_HM1 = cell(numPanels,1);
Rear_Irradiance_HM1_abs = cell(numPanels,1);

for k=1:numPanels
% --- Front Total
    Front_Irradiance{k} = IB_Front{k} + albedo_spec_per_PV{k} + ID_Iso_Front + ID_Cir_Front + ID_Hor_Front;
    Front_Irradiance_abs{k} = IB_Front_abs{k} + albedo_spec_per_PV_abs{k} + ID_Iso_Front_abs + ...
        ID_Cir_Front_abs + ID_Hor_Front_abs;

    Front_Irradiance_HM1{k} = IB_Front_HM1{k} + albedo_spec_per_HM1{k} + ID_Iso_Front_HM1 + ...
        ID_Cir_Front_HM1 + ID_Hor_Front_HM1;
    Front_Irradiance_HM1_abs{k} = IB_Front_abs_HM1{k} + albedo_spec_per_HM1_abs{k} + ID_Iso_Front_abs_HM1 + ...
        ID_Cir_Front_abs_HM1 + ID_Hor_Front_abs_HM1;

    Rear_Irradiance_HM1{k} = IB_Back_HM1{k} + albedo_spec_per_HM1_rear{k} +...
        ID_Iso_Rear_HM1 + ...
        ID_Cir_Rear_HM1 + ID_Hor_Rear_HM1;
    Rear_Irradiance_HM1_abs{k} = IB_Back_HM1_abs{k} +...
        albedo_spec_per_HM1_rear_abs{k} + ID_Iso_Rear_abs_HM1 + ...
        ID_Cir_Rear_abs_HM1 + ID_Hor_Rear_abs_HM1;
end

% --- Direct Beam (Rear)
% Preallocate cell array
IB_Back = cell(numPanelRows*numPanelColumns, 1);
IB_Back_abs = cell(numPanelRows*numPanelColumns, 1);

IB_Front_HM2 = cell(numPanelRows*numPanelColumns, 1);
IB_Front_abs_HM2 = cell(numPanelRows*numPanelColumns, 1);
IB_Back_HM2 = cell(numPanelRows*numPanelColumns, 1);
IB_Back_HM2_abs = cell(numPanelRows*numPanelColumns, 1);


for p = 1:numPanelRows*numPanelColumns
    
     IB_Back{p} = DNI.* cosd(AOI_HM1.Back) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  R_HM1.Back.Beam.*dirHorznIrradNormHourlyMod;
     IB_Back_abs{p} = DNI.* cosd(AOI_HM1.Back) .* ...
                  (1 - shadingFractionsHM1_direct(p)) .* ...
                  A_Panel.Back.Beam.* ...
                  R_HM1.Back.Beam.*dirHorznIrradNormHourlyMod;

     IB_Front_HM2{p} = DNI.* cosd(AOI_HM2.Front) .* ...
                  (1 - shadingFractionsHM2_direct(p)) .* ...
                  dirHorznIrradNormHourlyMod;
     IB_Front_abs_HM2{p} = DNI.* cosd(AOI_HM2.Front) .* ...
                  (1 - shadingFractionsHM2_direct(p)) .* ...
                  A_HM2.Front.Beam.* ...
                  dirHorznIrradNormHourlyMod;

     IB_Back_HM2{p} = DNI.* cosd(AOI_HM2.Back) .* ...
                  (1 - shadingFractionsHM2_direct(p)) .* ...
                  dirHorznIrradNormHourlyMod;
     IB_Back_HM2_abs{p} = DNI.* cosd(AOI_HM2.Back) .* ...
                  (1 - shadingFractionsHM2_direct(p)) .* ...
                  A_HM2.Back.Beam.* ...
                 dirHorznIrradNormHourlyMod;
   
     IB_Back{p}(IB_Back{p} < 0) = 0;    
     IB_Back_abs{p}(IB_Back_abs{p} < 0) = 0;     

     IB_Back_HM2{p}(IB_Back_HM2{p} < 0) = 0;    
     IB_Back_HM2_abs{p}(IB_Back_HM2_abs{p} < 0) = 0;

end



% --- Sky Diffuse (Rear)
[~, Iso_Front, Cir_Front, Hor_Front] = pvl_perez(hm2Tilt, hm2Azimuth, ...
    DHI, DNI, HExtra, zenith, azimuth, AM,difHorznIrradNormHourlyMod,dirHorznIrradNormHourlyMod,wavelengthDir, '1990');
[~, Iso_Rear, Cir_Rear, Hor_Rear] = pvl_perez(180-hm2Tilt, ...
    mod(hm2Azimuth + 180, 360), ...
    DHI, DNI, HExtra, zenith, azimuth, AM,difHorznIrradNormHourlyMod,...
    dirHorznIrradNormHourlyMod,wavelengthDir, '1990');

VFHMPanel = VF_HM_to_Panel{1,2};
ID_Iso_Rear = Iso_Front .* ...
               R_HM2.Front.Iso .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Cir_Rear = Cir_Front .* ...
               R_HM2.Front.Beam .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Hor_Rear = Hor_Front .* ...
               R_HM2.Front.Hor .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Iso_Rear_abs = Iso_Front .* ...
               R_HM2.Front.Iso .* ...
               reduction_dualHM2 .*A_Panel.Back.Iso.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);
ID_Cir_Rear_abs = Cir_Front .* ...
               R_HM2.Front.Beam .* ...
               reduction_dualHM2 .*A_Panel.Back.Beam.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);
ID_Hor_Rear_abs = Hor_Front .* ...
               R_HM2.Front.Hor .* ...
               reduction_dualHM2 .*A_Panel.Back.Hor.* ...
               difHorznIrradNormHourlyMod*max([VFHMPanel{:}]);

ID_Iso_Front_HM2 = Iso_Front .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Cir_Front_HM2 = Cir_Front .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Hor_Front_HM2 = Hor_Front .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Iso_Front_abs_HM2 = Iso_Front .* ...
               reduction_dualHM2 .*A_HM2.Front.Iso.* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Front_abs_HM2 = Cir_Front .* ...
               reduction_dualHM2 .*A_HM2.Front.Beam.* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Front_abs_HM2 = Hor_Front .* ...
               reduction_dualHM2 .*A_HM2.Front.Hor.* ...
               difHorznIrradNormHourlyMod;

ID_Iso_Rear_HM2 = Iso_Rear .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Cir_Rear_HM2 = Cir_Rear .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Hor_Rear_HM2 = Hor_Rear .* ...
               reduction_dualHM2 .* ...
               difHorznIrradNormHourlyMod;

ID_Iso_Rear_abs_HM2 = Iso_Rear .* ...
               reduction_dualHM2 .*A_HM2.Back.Iso.* ...
               difHorznIrradNormHourlyMod;
ID_Cir_Rear_abs_HM2 = Cir_Rear .* ...
               reduction_dualHM2 .*A_HM2.Back.Beam.* ...
               difHorznIrradNormHourlyMod;
ID_Hor_Rear_abs_HM2 = Hor_Rear .* ...
               reduction_dualHM2 .*A_HM2.Back.Hor.* ...
               difHorznIrradNormHourlyMod;


% % --- Rear Albedo
% I_Alb_Rear = pvl_Purdue_albedo_model(panelTilt_back, panelAzimuth_back, ...
%     panelElevation / panelHeight, groundAlbedo, ...
%     DHI, DNI, HExtra, zenith, azimuth, AM, ...
%     shadingFractions_direct, reduction_dif, VF_all, VF_ground,albedo_spec,...
%     albedo_spec,wavelengthDir, '1990');
I_Alb_Rear = albedo_per_PV_rear;
I_Alb_Rear(I_Alb_Rear < 0) = 0;

I_Alb_Rear_abs = albedo_per_PV_rear_abs;
I_Alb_Rear_abs(I_Alb_Rear_abs < 0) = 0;

I_Alb_HM2_Front = albedo_per_HM2;
I_Alb_HM2_Front(I_Alb_HM2_Front < 0) = 0;
I_Alb_HM2_Front_abs = albedo_per_HM2_abs;
I_Alb_HM2_Front_abs(I_Alb_HM2_Front_abs < 0) = 0;

I_Alb_HM2_Rear = albedo_per_HM2_rear;
I_Alb_HM2_Rear(I_Alb_HM2_Rear < 0) = 0;
I_Alb_HM2_Rear_abs = albedo_per_HM2_rear_abs;
I_Alb_HM2_Rear_abs(I_Alb_HM2_Rear_abs < 0) = 0;
%I_Alb_Rear = reshape(I_Alb_Rear, numPanelColumns, numPanelRows)';  % 4×3 matrix (rows x cols)
% --- Rear Total
Rear_Irradiance = cell(numPanels,1);
Rear_Irradiance_abs = cell(numPanels,1);

Front_Irradiance_HM2 = cell(numPanels,1);
Front_Irradiance_HM2_abs = cell(numPanels,1);

Rear_Irradiance_HM2 = cell(numPanels,1);
Rear_Irradiance_HM2_abs = cell(numPanels,1);

for k=1:numPanels
% --- Front Total
    Rear_Irradiance{k} = IB_Back{k} + albedo_spec_per_PV_rear{k} + ID_Iso_Rear + ID_Cir_Rear + ID_Hor_Rear;

    Rear_Irradiance_abs{k} = IB_Back_abs{k} + albedo_spec_per_PV_rear_abs{k} ...
        + ID_Iso_Rear_abs + ID_Cir_Rear_abs + ID_Hor_Rear_abs;

    Front_Irradiance_HM2{k} = IB_Front_HM2{k} + albedo_spec_per_HM2{k} + ...
        ID_Iso_Front_HM2 + ID_Cir_Front_HM2 + ID_Hor_Front_HM2;

    Front_Irradiance_HM2_abs{k} = IB_Front_abs_HM2{k} + albedo_spec_per_HM2_abs{k} ...
        + ID_Iso_Front_abs_HM2 + ID_Cir_Front_abs_HM2 + ID_Hor_Front_abs_HM2;

     Rear_Irradiance_HM2{k} = IB_Back_HM2{k} + albedo_spec_per_HM2_rear{k} + ...
        ID_Iso_Rear_HM2 + ID_Cir_Rear_HM2 + ID_Hor_Rear_HM2;

    Rear_Irradiance_HM2_abs{k} = IB_Back_HM2_abs{k} + albedo_spec_per_HM2_abs{k} ...
        + ID_Iso_Rear_abs_HM2 + ID_Cir_Rear_abs_HM2 + ID_Hor_Rear_abs_HM2;
end

% --- Output structure
out = struct();

% % View factors
% out.VF_all   = VF_all;
% out.VF_under = VF_under;
% out.VF_alley = VF_alley;
% out.VF_ground = VF_ground;

% Front components
out.IB_Front       = IB_Front;
out.IB_Front_HM1 = IB_Front_HM1;
out.IB_Front_HM2 = IB_Front_HM2;
out.I_Alb_Front    = I_Alb_Front;
out.I_Alb_HM1_Front = I_Alb_HM1_Front;
out.I_Alb_HM2_Front = I_Alb_HM2_Front;
out.ID_Iso_Front   = ID_Iso_Front;
out.ID_Iso_Front_HM1 = ID_Iso_Front_HM1;
out.ID_Iso_Front_HM2 = ID_Iso_Front_HM2;
out.ID_Cir_Front   = ID_Cir_Front;
out.ID_Cir_Front_HM1 = ID_Cir_Front_HM1;
out.ID_Cir_Front_HM2 = ID_Cir_Front_HM2;
out.ID_Hor_Front   = ID_Hor_Front;
out.ID_Hor_Front_HM1 = ID_Hor_Front_HM1;
out.ID_Hor_Front_HM2 = ID_Hor_Front_HM2;
out.IB_Front_abs=IB_Front_abs;
out.IB_Front_abs_HM1 = IB_Front_abs_HM1;
out.IB_Front_abs_HM2 = IB_Front_abs_HM2;
out.ID_Iso_Front_abs=ID_Iso_Front_abs;
out.ID_Iso_Front_abs_HM1 = ID_Iso_Front_abs_HM1;
out.ID_Iso_Front_abs_HM2 = ID_Iso_Front_abs_HM2;
out.ID_Cir_Front_abs=ID_Cir_Front_abs;
out.ID_Cir_Front_abs_HM1 = ID_Cir_Front_abs_HM1;
out.ID_Cir_Front_abs_HM2 = ID_Cir_Front_abs_HM2;
out.ID_Hor_Front_abs=ID_Hor_Front_abs;
out.ID_Hor_Front_abs_HM1 = ID_Hor_Front_abs_HM1;
out.ID_Hor_Front_abs_HM2 = ID_Hor_Front_abs_HM2;
out.albedo_spec_per_PV_abs=albedo_spec_per_PV_abs;
out.albedo_spec_per_HM1_abs=albedo_spec_per_HM1_abs;
out.albedo_spec_per_HM2_abs=albedo_spec_per_HM2_abs;
out.albedo_per_PV_abs=albedo_per_PV_abs;
out.albedo_per_HM1_Front_abs=albedo_per_HM1_abs;
out.albedo_per_HM2_Front_abs=albedo_per_HM2_abs;
out.I_Alb_Front_abs=I_Alb_Front_abs;
out.I_Alb_HM1_Front_abs=I_Alb_HM1_Front_abs;
out.I_Alb_HM2_Front_abs=I_Alb_HM2_Front_abs;
out.Front_Irradiance_abs=Front_Irradiance_abs;
out.Front_Irradiance_HM1_abs=Front_Irradiance_HM1_abs;
out.Front_Irradiance_HM2_abs=Front_Irradiance_HM2_abs;

% Rear components
out.IB_Back       = IB_Back;
out.IB_Back_HM1 = IB_Back_HM1;
out.IB_Back_HM2 = IB_Back_HM2;
out.I_Alb_Rear    = I_Alb_Rear;
out.I_Alb_HM1_Rear = I_Alb_HM1_Rear;
out.I_Alb_HM2_Rear = I_Alb_HM2_Rear;
out.ID_Iso_Rear   = ID_Iso_Rear;
out.ID_Iso_Rear_HM1 = ID_Iso_Rear_HM1;
out.ID_Iso_Rear_HM2 = ID_Iso_Rear_HM2;
out.ID_Cir_Rear   = ID_Cir_Rear;
out.ID_Cir_Rear_HM1=ID_Cir_Rear_HM1;
out.ID_Cir_Rear_HM2=ID_Cir_Rear_HM2;
out.ID_Hor_Rear   = ID_Hor_Rear;
out.ID_Hor_Rear_HM1=ID_Hor_Rear_HM1;
out.ID_Hor_Rear_HM2=ID_Hor_Rear_HM2;
out.IB_Back_abs=IB_Back_abs;
out.IB_Back_HM1_abs=IB_Back_HM1_abs;
out.IB_Back_HM2_abs=IB_Back_HM2_abs;
out.ID_Iso_Rear_abs=ID_Iso_Rear_abs;
out.ID_Iso_Rear_abs_HM1=ID_Iso_Rear_abs_HM1;
out.ID_Iso_Rear_abs_HM2=ID_Iso_Rear_abs_HM2;
out.ID_Cir_Rear_abs=ID_Cir_Rear_abs;
out.ID_Cir_Rear_abs_HM1=ID_Cir_Rear_abs_HM1;
out.ID_Cir_Rear_abs_HM2=ID_Cir_Rear_abs_HM2;
out.ID_Hor_Rear_abs=ID_Hor_Rear_abs;
out.ID_Hor_Rear_abs_HM1=ID_Hor_Rear_abs_HM1;
out.ID_Hor_Rear_abs_HM2=ID_Hor_Rear_abs_HM2;
out.I_Alb_Rear_abs=I_Alb_Rear_abs;
out.I_Alb_HM1_Rear_abs=I_Alb_HM1_Rear_abs;
out.I_Alb_HM2_Rear_abs=I_Alb_HM2_Rear_abs;
out.albedo_spec_per_PV_rear_abs=albedo_spec_per_PV_rear_abs;
out.albedo_spec_per_HM1_rear_abs=albedo_spec_per_HM1_rear_abs;
out.albedo_spec_per_HM2_rear_abs=albedo_spec_per_HM2_rear_abs;
out.albedo_per_PV_rear_abs=albedo_per_PV_rear_abs;
out.albedo_per_HM1_rear_abs=albedo_per_HM1_rear_abs;
out.albedo_per_HM2_rear_abs=albedo_per_HM2_rear_abs;
out.Back_Irradiance_abs=Rear_Irradiance_abs;
out.Rear_Irradiance_HM1_abs=Rear_Irradiance_HM1_abs;
out.Rear_Irradiance_HM2_abs=Rear_Irradiance_HM2_abs;


% --- Clean non-numeric outputs ---
outFields = fieldnames(out);

for f = 1:numel(outFields)
    val = out.(outFields{f});

    % Case 1: if it's a cell array, ensure each cell is numeric
    if iscell(val)
        for k = 1:numel(val)
            if ~isnumeric(val{k}) || any(isnan(val{k}(:))) || isempty(val{k})
                valSize = size(val{k});
                if isempty(valSize)
                    valSize = [1 1];
                end
                val{k} = zeros(valSize);
            end
        end
        out.(outFields{f}) = val;

    % Case 2: numeric array — set NaNs to 0
    elseif isnumeric(val)
        val(isnan(val)) = 0;
        out.(outFields{f}) = val;

    % Case 3: anything else (string, struct, etc) — replace with 0
    else
        %warning(['Non-numeric field "', outFields{f}, '" replaced with 0.']);
        out.(outFields{f}) = 0;
    end
end


end
