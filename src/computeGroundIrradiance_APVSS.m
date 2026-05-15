function [PAR_inten, ePAR_inten,IR_inten, Ri_inten,PAR_inten_mmole, ...
    ePAR_inten_mmole,Ri_spec, shadingMaskhm1,...
    shadingMaskhm2,directSolar,diffuseSolar,idx_PAR,idx_ePAR] = ...
    computeGroundIrradiance_APVSS(...
    groundX, groundY, DNI,IB_Front_HM1, IB_Front_HM2, ID_Iso_Front_HM1, ID_Iso_Front_HM2,...
        ID_Cir_Front_HM1, ID_Cir_Front_HM2,...
        ID_Hor_Front_HM1, ID_Hor_Front_HM2, ...
        theta, shadowPolygons, ...
        viewFactorCell_GroundToHM1Rear, viewFactorCell_GroundToHM2Rear, ...
        panel_reflected_front, panel_reflected_rear, ...
        dirHorznIrradNormHourlyMod, difHorznIrradNormHourlyMod, wavelengthDir,...
        T_HM1,T_HM2,ID_Iso_sky, ID_Cir_sky, ID_Hor_sky,SVF_ground)



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

    %% --- Compute Shading Mask ---
    shadingMaskhm1 = zeros(size(groundX));  % 0 = lit; 1 = shaded
    shadingMaskhm2 = zeros(size(groundX));
    %% --- Preallocate Matrices for Integrated Intensities ---
    PAR_inten_mmole  = zeros(numGroundY, numGroundX);  % in μmol m⁻² s⁻¹
    ePAR_inten_mmole = zeros(numGroundY, numGroundX);  % in μmol m⁻² s⁻¹
    PAR_inten  = zeros(numGroundY, numGroundX);  % in W m-2
    ePAR_inten = zeros(numGroundY, numGroundX);  % in W m-2
    IR_inten   = zeros(numGroundY, numGroundX);  % in W m-2
    Ri_inten   = zeros(numGroundY, numGroundX);  % in W m-2 
    directSolar = cell(numGroundY, numGroundX);  % in W m-2 
    diffuseSolar = cell(numGroundY, numGroundX); % in W m-2 
    Ri_spec = cell(numGroundY, numGroundX); % in W m-2 
    %% --- Define Wavelength Indices for Each Band ---
         idx_PAR  = find(wavelengthDir >= 400 & wavelengthDir <= 700);
         idx_ePAR = find(wavelengthDir >= 400 & wavelengthDir <= 750);
         idx_IR   = find(wavelengthDir >= 750 & wavelengthDir <= 4000);
         idx_ri = find(wavelengthDir >= 285 & wavelengthDir <= 800);
     numPanels = length(shadowPolygons);
     direct_component = cell(numGroundY, numGroundX);
     diffuse_component = cell(numGroundY, numGroundX);
    for idx = 1:numPanels
        polyXhm1 = shadowPolygons{idx}{1,1}(:,1);
        polyYhm1 = shadowPolygons{idx}{1,1}(:,2);
        inPolyhm1 = inpolygon(groundX(:), groundY(:), polyXhm1, polyYhm1);
        shadingMaskhm1(inPolyhm1) = 1;
        polyXhm2 = shadowPolygons{idx}{1,2}(:,1);
        polyYhm2 = shadowPolygons{idx}{1,2}(:,2);
        inPolyhm2 = inpolygon(groundX(:), groundY(:), polyXhm2, polyYhm2);
        shadingMaskhm2(inPolyhm2) = 1;
        for i = 1:numGroundY
            for j = 1:numGroundX

                %VFf = VF_front_combined(i,j);
                %VFr = VF_rear_combined(i,j);
                % Determine spectral components based on shading
                if shadingMaskhm1(i,j) == 1
                    % Shaded — use panel transmission
                    if(trapz(wavelengthDir,IB_Front_HM1{idx})>0)
                        direct_component{i,j} = IB_Front_HM1{idx}.* T_HM1.Front.Beam;
                    else
                        direct_component{i,j}=zeros(size(wavelengthDir,1),1);
                    end
                        
                else
                    if shadingMaskhm2(i,j) == 1
                        if(trapz(wavelengthDir,IB_Front_HM2{idx})>0)
                            direct_component{i,j} = IB_Front_HM2{idx}.* T_HM2.Front.Beam;
                        else
                            direct_component{i,j}=zeros(size(wavelengthDir,1),1);
                        end

                                                                         
                    else


                    % Unshaded — full sky view
                        direct_component{i,j} = Idirect_ground;
                   
                    end


                end

                 
            end

        end

    end

    for i=1:numGroundY
        for j=1:numGroundX
            diffuse_component{i,j} = (((ID_Iso_Front_HM1.*T_HM1.Front.Iso)+...
                        (ID_Cir_Front_HM1.*T_HM1.Front.Beam) + ...
                        (ID_Hor_Front_HM1.*T_HM1.Front.Hor)).* ...
                        (1-SVF_ground(i,j)))+ ...
                        (((ID_Iso_sky + ID_Cir_sky + ID_Hor_sky) .* SVF_ground(i,j)).*difHorznIrradNormHourlyMod);
            
                 % Always include reflection
                reflected_component = panel_reflected_front{i,j} + panel_reflected_rear{i,j};

                % Combine all sources
                spec = direct_component{i,j} + diffuse_component{i,j} + reflected_component;
                directSolar{i,j} = direct_component{i,j}; %Need to multiply with extinction coefficient to assess plant factors
                diffuseSolar{i,j} = diffuse_component{i,j} + reflected_component;
                Ri_spec{i,j} = spec(idx_ri);

                %--- Integrate Over the Bands ---
                % For PAR and ePAR, convert integrated energy (W/m²) into photon flux (μmol m⁻² s⁻¹)
                % using:   photon_flux = (1e6/NA) * ∫[ I(λ)*(λ*1e-9)/(h*c) dλ ]
                %
                % For Infrared, simply integrate the spectral irradiance over 750–4000 nm.
            
                PAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) .* ((wavelengthDir(idx_PAR)*1e-9)/(h*c)) );
                          
                ePAR_inten_mmole(i,j) = (1e6/NA) * trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) .* ((wavelengthDir(idx_ePAR)*1e-9)/(h*c)) );
                          
                IR_inten(i,j)  = trapz( wavelengthDir(idx_IR), spec(idx_IR) );

                Ri_inten(i,j)  = trapz( wavelengthDir(idx_ri), spec(idx_ri) );

                PAR_inten(i,j) =  trapz( wavelengthDir(idx_PAR), ...
                              spec(idx_PAR) );
                          
                ePAR_inten(i,j) = trapz( wavelengthDir(idx_ePAR), ...
                              spec(idx_ePAR) );
        end
    end
    a=1;

    %end
    % Optionally, you might display a note:
    % disp('Integrated intensities (PAR, ePAR, IR) computed for each ground cell.');
end
