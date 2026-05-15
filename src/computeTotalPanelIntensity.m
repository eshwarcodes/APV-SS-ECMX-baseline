function [totalFrontIrradiance, totalRearIrradiance, reflectedFrontIntensity, reflectedRearIntensity] = ...
    computeTotalPanelIntensity(Front_Irradiance, Rear_Irradiance, panel_reflected_front, panel_reflected_rear, wavelength)
% computeTotalIntensity computes the total integrated intensity from spectral
% measurements.
%
% The function integrates over the spectral range (given by the wavelength vector)
% for the following:
%   1. Front_Irradiance: cell array (numPanels x 1) each containing an 1882x1 spectrum.
%   2. Rear_Irradiance: cell array (numPanels x 1) each containing an 1882x1 spectrum.
%   3. panel_reflected_front: cell array (e.g., 200x200) each containing an 1882x1 spectrum.
%   4. panel_reflected_rear: cell array (e.g., 200x200) each containing an 1882x1 spectrum.
%
% The integration is done using the trapezoidal rule (trapz).
%
% Inputs:
%   Front_Irradiance         - Cell array (numPanels x 1), each cell an 1882x1 vector (W/m²/nm)
%   Rear_Irradiance          - Cell array (numPanels x 1), each cell an 1882x1 vector (W/m²/nm)
%   panel_reflected_front    - Cell array (ground grid size, e.g. 200x200), each cell an 1882x1 vector (W/m²/nm)
%   panel_reflected_rear     - Cell array (ground grid size, e.g. 200x200), each cell an 1882x1 vector (W/m²/nm)
%   wavelength               - Vector (1882x1) of wavelengths (nm) corresponding to the spectral data
%
% Outputs:
%   totalFrontIrradiance         - Vector (numPanels x 1) integrated intensity from Front_Irradiance (W/m²)
%   totalRearIrradiance          - Vector (numPanels x 1) integrated intensity from Rear_Irradiance (W/m²)
%   totalPanelReflectedIrradiance- Matrix (same size as panel_reflected_front cell array) of the sum of
%                                  integrated reflected intensity from the front and rear sides (W/m²)
%

    %% --- Integrate Front and Rear Irradiance for Panels ---
    % Use cellfun to apply trapz to each cell.
    totalFrontIrradiance = cellfun(@(spec) trapz(double(wavelength), double(spec)), Front_Irradiance);
    totalRearIrradiance  = cellfun(@(spec) trapz(double(wavelength), double(spec)), Rear_Irradiance);
    
    %% --- Integrate Reflected Irradiance on the Ground ---
    % For the ground, each cell in panel_reflected_front and panel_reflected_rear is integrated.
    reflectedFrontIntensity = cellfun( ...
    @(spec) trapz(double(wavelength), double(spec)), ...
    panel_reflected_front);
    reflectedRearIntensity = cellfun( ...
    @(spec) trapz(double(wavelength), double(spec)), ...
    panel_reflected_rear);
    
  
    
end
