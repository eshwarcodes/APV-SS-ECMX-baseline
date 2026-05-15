function [sun_interp,shade_interp,...
    sun_albedo, shade_albedo, ...
    dirHorznIrradNormHourlyMod_new, difHorznIrradNormHourlyMod_new] = ...
    calculate_spectral_components(dirHorznIrradNormHourlyMod, difHorznIrradNormHourlyMod, ...
    wavelengthDir, wavelength,sun_tables,shade_tables)


% Integrate spectra
integdirnorm_temp = trapz(wavelengthDir, dirHorznIrradNormHourlyMod);
integdifnorm_temp = trapz(wavelengthDir,    difHorznIrradNormHourlyMod);

% Safe normalization: if integral = 0, result = all zeros
if integdirnorm_temp == 0
    dirHorznIrradNormHourlyMod_new = zeros(size(dirHorznIrradNormHourlyMod));
else
    dirHorznIrradNormHourlyMod_new = dirHorznIrradNormHourlyMod / integdirnorm_temp;
end

if integdifnorm_temp == 0
    difHorznIrradNormHourlyMod_new = zeros(size(difHorznIrradNormHourlyMod));
else
    difHorznIrradNormHourlyMod_new = difHorznIrradNormHourlyMod / integdifnorm_temp;
end


% Initialize output structs
sun_interp = struct();
shade_interp = struct();
sun_albedo = struct();
shade_albedo = struct();

% Loop over LAI levels
for i = 1:length(sun_tables)
    % ---------- SUNLIT ----------
    T_sun = sun_tables{i};

    sun_interp(i).CanopyReflectance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.CanopyReflectance, wavelengthDir, 'linear', 'extrap'));
    sun_interp(i).CanopyTransmittance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.CanopyTransmittance, wavelengthDir, 'linear', 'extrap'));
    sun_interp(i).CanopyAbsorptance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.CanopyAbsorptance, wavelengthDir, 'linear', 'extrap'));
    sun_interp(i).LeafReflectance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.LeafReflectance, wavelengthDir, 'linear', 'extrap'));
    sun_interp(i).LeafTransmittance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.LeafTransmittance, wavelengthDir, 'linear', 'extrap'));
    sun_interp(i).LeafAbsorptance = clamp01(interp1(T_sun.Wavelength_nm, T_sun.LeafAbsorptance, wavelengthDir, 'linear', 'extrap'));

    % Compute normalized spectral albedo
    spec = sun_interp(i).CanopyReflectance .* (dirHorznIrradNormHourlyMod + difHorznIrradNormHourlyMod);
    sun_albedo(i).Spectral = spec ./ trapz(wavelengthDir, spec);

    % ---------- SHADED ----------
    T_shade = shade_tables{i};

    shade_interp(i).CanopyReflectance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.CanopyReflectance, wavelengthDir, 'linear', 'extrap'));
    shade_interp(i).CanopyTransmittance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.CanopyTransmittance, wavelengthDir, 'linear', 'extrap'));
    shade_interp(i).CanopyAbsorptance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.CanopyAbsorptance, wavelengthDir, 'linear', 'extrap'));
    shade_interp(i).LeafReflectance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.LeafReflectance, wavelengthDir, 'linear', 'extrap'));
    shade_interp(i).LeafTransmittance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.LeafTransmittance, wavelengthDir, 'linear', 'extrap'));
    shade_interp(i).LeafAbsorptance = clamp01(interp1(T_shade.Wavelength_nm, T_shade.LeafAbsorptance, wavelengthDir, 'linear', 'extrap'));

    % Compute normalized spectral albedo
    spec = shade_interp(i).CanopyReflectance .* (dirHorznIrradNormHourlyMod + difHorznIrradNormHourlyMod);
    shade_albedo(i).Spectral = spec ./ trapz(wavelengthDir, spec);
end



end

% --- Helper function ---
function val = get_stack_val(AOI, stack_allangles, angle_grid)
    if AOI == 90
        val = [1 0 0]; % [R T A] dummy values
        val = val(1) * ones(1, 341);
    else
        val = interp1(angle_grid, stack_allangles, AOI, 'linear');
    end
end
