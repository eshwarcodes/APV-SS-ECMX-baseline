function [PAR_inten, ePAR_inten, IR_inten,PAR_inten_mmole, ePAR_inten_mmole, IR_inten_mmole, shadingMask,...
    directSolar,diffuseSolar,idx_PAR,idx_ePAR] = ...
    computeGroundIrradiance(...
    groundX, groundY, DNI,ID_Iso_Front,ID_Cir_Front,ID_Hor_Front,theta, shadowPolygons, ...
    viewFactor_RearPVGround, viewFactorCell_GroundToPVRear, ...
    panel_reflected_front, panel_reflected_rear, ...
    dirHorznIrradNormHourlyMod, difHorznIrradNormHourlyMod, wavelengthDir,IB_Front,...
    T_Beam_front_interp ,T_Iso_front_interp ,T_Hor_front_interp,ID_Iso_sky, ID_Cir_sky, ID_Hor_sky)
% computeGroundIrradiance computes spatial ground irradiance over a spectral
% range (280–4000 nm) and then integrates it into three bands:
%   PAR (400–700 nm), ePAR (400–750 nm) and Infrared (750–4000 nm).
%
% The PAR and ePAR values (integrated from spectral irradiance in W/m2) are
% converted into photon flux (μmol m⁻² s⁻¹) using the conversion:
%
%   Φ = (1e6/NA) * ∫[I(λ) * ((λ*1e-9)/(h*c)) dλ]
%
% where h is Planck's constant, c is the speed of light and NA is Avogadro's 
% constant.
%
% Inputs:
%   groundX, groundY          - Meshgrid coordinates of the ground (e.g. 200x200)
%   DNI, DHI                  - Direct Normal and Diffuse Horizontal Irradiance
%   theta                     - Solar zenith angle (radians)
%   shadowPolygons            - Cell array of vertices for polygons casting shadows
%   viewFactor_FrontPVGround  - Cell array of sky view factor maps affected by PV front
%   viewFactor_RearPVGround   - Cell array of sky view factor maps affected by PV rear
%   panel_reflected_front     - Cell array of spectral irradiance from front side reflections
%   panel_reflected_rear      - Cell array of spectral irradiance from rear side reflections
%   dirHorznIrradNormHourlyMod - Modifier for direct irradiance on horizontal ground
%   difHorznIrradNormHourlyMod  - Modifier for diffuse irradiance on horizontal ground
%   wavelengthDir             - Vector of wavelengths (nm), e.g. 280:5:4000
%
% Outputs:
%   PAR_inten  - Matrix (numGroundY x numGroundX) of integrated PAR intensity 
%                in μmol m⁻² s⁻¹ (400–700 nm)
%   ePAR_inten - Matrix (numGroundY x numGroundX) of integrated extended PAR intensity 
%                in μmol m⁻² s⁻¹ (400–750 nm)
%   IR_inten   - Matrix (numGroundY x numGroundX) of integrated Infrared irradiance (750–4000 nm) in W/m²
%   shadingMask- Binary shading mask (1: shaded, 0: lit)
%

    % Determine grid size
    [numGroundY, numGroundX] = size(groundX);
    
    %% Constants for photon conversion
    h  = 6.62607015e-34;    % Planck's constant (J·s)
    c  = 2.99792458e8;      % Speed of light (m/s)
    NA = 6.02214076e23;     % Avogadro's number (1/mol)

    %% --- Direct and Diffuse Irradiance on Horizontal Ground ---
    Idirect_ground = DNI * cos(theta) * dirHorznIrradNormHourlyMod;  % in W/m²
    %Idiff_ground   = DHI * difHorznIrradNormHourlyMod;               % in W/m²/nm

    %% --- Compute Shading Mask ---
    shadingMask = zeros(size(groundX));  % 0 = lit; 1 = shaded
    %% --- Preallocate Matrices for Integrated Intensities ---
    PAR_inten_mmole  = zeros(numGroundY, numGroundX);  % in μmol m⁻² s⁻¹
    ePAR_inten_mmole = zeros(numGroundY, numGroundX);  % in μmol m⁻² s⁻¹
    IR_inten_mmole   = zeros(numGroundY, numGroundX);  % in μmol m⁻² s⁻¹

    PAR_inten  = zeros(numGroundY, numGroundX);  % in W m-2
    ePAR_inten = zeros(numGroundY, numGroundX);  % in W m-2
    IR_inten   = zeros(numGroundY, numGroundX);  % in W m-2
    directSolar = cell(numGroundY, numGroundX);
    diffuseSolar = cell(numGroundY, numGroundX);
    
    %% --- Define Wavelength Indices for Each Band ---
    idx_PAR  = find(wavelengthDir >= 400 & wavelengthDir <= 700);
    idx_ePAR = find(wavelengthDir >= 400 & wavelengthDir <= 750);
    idx_IR   = find(wavelengthDir >= 750 & wavelengthDir <= 4000);
    numPanels = length(shadowPolygons); 
    for idx = 1:numPanels
        polyX = shadowPolygons{idx}(:,1);
        polyY = shadowPolygons{idx}(:,2);
        inPoly = inpolygon(groundX(:), groundY(:), polyX, polyY);
        shadingMask(inPoly) = 1;
    end

    
    %% --- Loop over Each Ground Cell, Sum Spectral Contributions and Integrate ---
    % The spectral irradiance for a cell is computed as the sum of:
    % 1. Direct irradiance (only if not shaded)
    % 2. Diffuse effective irradiance from the sky
    % 3. Reflected irradiance from the PV panels (front and rear)
    %
    % Note: panel_reflected_front and panel_reflected_rear are assumed to be cell arrays
    %       with each cell containing a spectrum (vector of length equal to wavelengthDir).
    
        for i = 1:numGroundY
            for j = 1:numGroundX
            
                if shadingMask(i,j) == 1
                    % Shaded — use panel transmission
                    direct_component = IB_Front{idx}.* T_Beam_front_interp;                                               
                else
                    % Unshaded — full sky view
                    direct_component = Idirect_ground;
                    
                end

            end
        end
    %end
    %% --- Combine Sky View Factors from All Panels ---
    %VF_front_combined = zeros(size(groundX));
    VF_rear_combined  = zeros(size(groundX));
    numPanels = length(viewFactor_RearPVGround);
    
    for idx = 1:numPanels
        %VF_front_combined = VF_front_combined + viewFactor_FrontPVGround{idx};
        VF_rear_combined  = VF_rear_combined  + viewFactorCell_GroundToPVRear{idx};
    end
    %VF_front_combined = VF_front_combined / numPanels;
    VF_rear_combined  = VF_rear_combined  / numPanels;

    VF_total = VF_rear_combined; 
    %+ ...
           % viewFactorCell_GroundToPVRear{1};

    % Normalize only where total > 1
    mask_over = VF_total > 1;
    VF_panel = VF_rear_combined;
    

    % Normalize the overshoot proportionally
    VF_panel(mask_over) = VF_panel(mask_over) ./ VF_total(mask_over);
    
    % Recompute total and SVF
    SVF_ground = 1 - (VF_panel);
    % + VF_PV
    SVF_ground = max(0, SVF_ground);  % optional: enforce floor

    for i = 1:numGroundY
            for j = 1:numGroundX  
                if shadingMask(i,j) == 1
                    diffuse_component = (((ID_Iso_Front.*T_Iso_front_interp)+...
                        (ID_Cir_Front.*T_Beam_front_interp) + ...
                        (ID_Hor_Front.*T_Hor_front_interp)).* ...
                        VF_panel(i,j))+ ...
                        (ID_Iso_sky + ID_Cir_sky + ID_Hor_sky) .* SVF_ground(i,j).*difHorznIrradNormHourlyMod;
                else
                     diffuse_component=(ID_Iso_sky + ID_Cir_sky + ID_Hor_sky).*difHorznIrradNormHourlyMod;
                end

                % Always include reflection
                reflected_component = panel_reflected_front{i,j} + panel_reflected_rear{i,j};

                % Combine all sources
                spec = direct_component + diffuse_component + reflected_component;
                directSolar{i}{j} = direct_component; %Need to multiply with extinction coefficient to assess plant factors
                diffuseSolar{i}{j} = diffuse_component + reflected_component;
            
                %--- Integrate Over the Bands ---
                % For PAR and ePAR, convert integrated energy (W/m²) into photon flux (μmol m⁻² s⁻¹)
                % using:   photon_flux = (1e6/NA) * ∫[ I(λ)*(λ*1e-9)/(h*c) dλ ]
                %
                % For Infrared, simply integrate the spectral irradiance over 750–4000 nm.
            
                PAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) .* ((wavelengthDir(idx_PAR)*1e-9)/(h*c)) );
                          
                ePAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) .* ((wavelengthDir(idx_ePAR)*1e-9)/(h*c)) );
                          
                IR_inten_mmole(i,j)  = trapz( wavelengthDir(idx_IR), spec(idx_IR) );

                PAR_inten(i,j) =  trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) );
                          
                ePAR_inten(i,j) = trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) );
                          
                IR_inten(i,j)  = trapz( wavelengthDir(idx_IR), spec(idx_IR) );

                %--- Integrate Over the Bands ---
                % For PAR and ePAR, convert integrated energy (W/m²) into photon flux (μmol m⁻² s⁻¹)
                % using:   photon_flux = (1e6/NA) * ∫[ I(λ)*(λ*1e-9)/(h*c) dλ ]
                %
                % For Infrared, simply integrate the spectral irradiance over 750–4000 nm.
            
                PAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) .* ((wavelengthDir(idx_PAR)*1e-9)/(h*c)) );
                          
                ePAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) .* ((wavelengthDir(idx_ePAR)*1e-9)/(h*c)) );
                          
                IR_inten_mmole(i,j)  = trapz( wavelengthDir(idx_IR), spec(idx_IR) );

                PAR_inten(i,j) =  trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) );
                          
                ePAR_inten(i,j) = trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) );
                          
                IR_inten(i,j)  = trapz( wavelengthDir(idx_IR), spec(idx_IR) );
            
            end

    end

       % a=1;
    % Optionally, you might display a note:
    % disp('Integrated intensities (PAR, ePAR, IR) computed for each ground cell.');
end
