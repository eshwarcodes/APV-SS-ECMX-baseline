function [AOI, AOI_Hor, R_interp, T_interp, A_interp] = computeHMOpticsDualFace( ...
    hmTilt, hmAzimuth, zenith, sunAz, ...
    R_stack_allangles, T_stack_allangles, ...
    A_stack_allangles, ...
    wavelength_optics, wavelengthDir)

    % FRONT FACE
    AOI.Front = pvl_getaoi(hmTilt, hmAzimuth, zenith, sunAz);
    AOI.Front((AOI.Front > 90) | (AOI.Front < 0)) = 90;

    AOI_Hor.Front = pvl_getaoi(hmTilt, hmAzimuth, 90, hmAzimuth);
    AOI_Hor.Front((AOI_Hor.Front > 90) | (AOI_Hor.Front < 0)) = 90;

    [R_Beam_F, R_Iso_F, R_Albedo_F] = ...
        pvl_iam_specular_spectral(hmTilt, AOI.Front, R_stack_allangles);
    [T_Beam_F, T_Iso_F, T_Albedo_F] = ...
        pvl_iam_specular_spectral_transmittance(hmTilt, AOI.Front, T_stack_allangles);
       
     [A_Beam_F, A_Iso_F, A_Albedo_F] = ...
        pvl_iam_specular_spectral_abs(hmTilt, AOI.Front, A_stack_allangles);

    R_Hor_F = pvl_iam_specular_spectral(hmTilt, AOI_Hor.Front, R_stack_allangles);
    T_Hor_F = pvl_iam_specular_spectral_transmittance(hmTilt, AOI_Hor.Front, T_stack_allangles); 
   
    A_Hor_F = pvl_iam_specular_spectral_abs(hmTilt, AOI_Hor.Front, A_stack_allangles);

    % BACK FACE
    hmTilt_back = 180 - hmTilt;
    hmAzimuth_back = mod(hmAzimuth + 180, 360);

    AOI.Back = pvl_getaoi(hmTilt_back, hmAzimuth_back, zenith, sunAz);
    AOI.Back((AOI.Back > 90) | (AOI.Back < 0)) = 90;

    AOI_Hor.Back = pvl_getaoi(hmTilt_back, hmAzimuth_back, 90, hmAzimuth_back);
    AOI_Hor.Back((AOI_Hor.Back > 90) | (AOI_Hor.Back < 0)) = 90;

    [R_Beam_B, R_Iso_B, R_Albedo_B] = ...
        pvl_iam_specular_spectral(hmTilt_back, AOI.Back, R_stack_allangles);
    [T_Beam_B, T_Iso_B, T_Albedo_B] = ...
    pvl_iam_specular_spectral_transmittance(hmTilt_back, AOI.Back, T_stack_allangles);     
    [A_Beam_B, A_Iso_B, A_Albedo_B] = ...
        pvl_iam_specular_spectral_abs(hmTilt_back, AOI.Back, A_stack_allangles);

    R_Hor_B = pvl_iam_specular_spectral(hmTilt_back, AOI_Hor.Back, R_stack_allangles);
    T_Hor_B = pvl_iam_specular_spectral_transmittance(hmTilt_back, AOI_Hor.Back, T_stack_allangles);
    
    A_Hor_B = pvl_iam_specular_spectral_abs(hmTilt_back, AOI_Hor.Back, A_stack_allangles);

    % Interpolate — Reflectance
    R_interp.Front.Beam   = clamp01(interp1(wavelength_optics, R_Beam_F, wavelengthDir, 'linear', 'extrap'));
    R_interp.Front.Iso    = clamp01(interp1(wavelength_optics, R_Iso_F, wavelengthDir, 'linear', 'extrap'));
    R_interp.Front.Albedo = clamp01(interp1(wavelength_optics, R_Albedo_F, wavelengthDir, 'linear', 'extrap'));
    R_interp.Front.Hor    = clamp01(interp1(wavelength_optics, R_Hor_F, wavelengthDir, 'linear', 'extrap'));

    R_interp.Back.Beam   = clamp01(interp1(wavelength_optics, R_Beam_B, wavelengthDir, 'linear', 'extrap'));
    R_interp.Back.Iso    = clamp01(interp1(wavelength_optics, R_Iso_B, wavelengthDir, 'linear', 'extrap'));
    R_interp.Back.Albedo = clamp01(interp1(wavelength_optics, R_Albedo_B, wavelengthDir, 'linear', 'extrap'));
    R_interp.Back.Hor    = clamp01(interp1(wavelength_optics, R_Hor_B, wavelengthDir, 'linear', 'extrap'));

    % Interpolate — Transmittance
    T_interp.Front.Beam   = clamp01(interp1(wavelength_optics, T_Beam_F, wavelengthDir, 'linear', 'extrap'));
    T_interp.Front.Iso    = clamp01(interp1(wavelength_optics, T_Iso_F, wavelengthDir, 'linear', 'extrap'));
    T_interp.Front.Albedo = clamp01(interp1(wavelength_optics, T_Albedo_F, wavelengthDir, 'linear', 'extrap'));
    T_interp.Front.Hor    = clamp01(interp1(wavelength_optics, T_Hor_F, wavelengthDir, 'linear', 'extrap'));

    T_interp.Back.Beam   = clamp01(interp1(wavelength_optics, T_Beam_B, wavelengthDir, 'linear', 'extrap'));
    T_interp.Back.Iso    = clamp01(interp1(wavelength_optics, T_Iso_B, wavelengthDir, 'linear', 'extrap'));
    T_interp.Back.Albedo = clamp01(interp1(wavelength_optics, T_Albedo_B, wavelengthDir, 'linear', 'extrap'));
    T_interp.Back.Hor    = clamp01(interp1(wavelength_optics, T_Hor_B, wavelengthDir, 'linear', 'extrap'));

    % Interpolate — Absorptance
    A_interp.Front.Beam   = clamp01(interp1(wavelength_optics, A_Beam_F, wavelengthDir, 'linear', 'extrap'));
    A_interp.Front.Iso    = clamp01(interp1(wavelength_optics, A_Iso_F, wavelengthDir, 'linear', 'extrap'));
    A_interp.Front.Albedo = clamp01(interp1(wavelength_optics, A_Albedo_F, wavelengthDir, 'linear', 'extrap'));
    A_interp.Front.Hor    = clamp01(interp1(wavelength_optics, A_Hor_F, wavelengthDir, 'linear', 'extrap'));

    A_interp.Back.Beam   = clamp01(interp1(wavelength_optics, A_Beam_B, wavelengthDir, 'linear', 'extrap'));
    A_interp.Back.Iso    = clamp01(interp1(wavelength_optics, A_Iso_B, wavelengthDir, 'linear', 'extrap'));
    A_interp.Back.Albedo = clamp01(interp1(wavelength_optics, A_Albedo_B, wavelengthDir, 'linear', 'extrap'));
    A_interp.Back.Hor    = clamp01(interp1(wavelength_optics, A_Hor_B, wavelengthDir, 'linear', 'extrap'));
end
