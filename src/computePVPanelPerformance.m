function PanelResult = computePVPanelPerformance( ...
    wavelengthDir, E0, ModuleParameters, ...
    wind_vel, Tao, numPanels, EgRef, dEgdt, ...
    Front_Irradiance_abs, bifaciality_factor, ...
    Rear_Irradiance_abs, bifacialModuleNamesSorted)

% Initialize output struct array
PanelResult(numPanels,1) = struct( ...
    'Tcell', [], 'Tmodule', [], ...
    'G_total', [], 'G_total_inten', [], ...
    'IL', [], 'I0', [], 'Rs', [], 'Rsh', [], 'nNsVth', [], ...
    'PCE', [], 'PCE_percent', [], ...
    'DiodeResult', []);  % 'DiodeResult' stores pvl_singlediode output

% Module area from database (assumes same for all panels)
A_module = bifacialModuleNamesSorted.A_c;

% Loop through each panel
for p = 1:numPanels
    % Combine front + bifacial rear spectrum
    G_total = Front_Irradiance_abs{p} + bifaciality_factor * Rear_Irradiance_abs{p};
    G_total_inten = trapz(wavelengthDir, G_total);  % integrated irradiance [W/m²]
    
    % Compute temperatures (SAPM)
    [Tcell, Tmodule] = pvl_sapmcelltemp( ...
        G_total_inten, E0, ...
        ModuleParameters.a, ModuleParameters.b, ...
        wind_vel, Tao, ModuleParameters.dT);

    % Compute diode model parameters
    [IL, I0, Rs, Rsh, nNsVth] = pvl_calcparams_desoto( ...
        G_total_inten, Tcell, ModuleParameters.alpha_isc, ...
        ModuleParameters, EgRef, dEgdt);

    % Solve IV curve using single diode model
    DiodeResult = pvl_singlediode(IL, I0, Rs, Rsh, nNsVth);

    % Compute PCE
    PCE = DiodeResult.Pmp / (G_total_inten * A_module);
    PCE_percent = 100 * PCE;

    % Store everything into output struct
    PanelResult(p).Tcell         = Tcell;
    PanelResult(p).Tmodule       = Tmodule;
    PanelResult(p).G_total       = G_total;
    PanelResult(p).G_total_inten = G_total_inten;
    PanelResult(p).IL            = IL;
    PanelResult(p).I0            = I0;
    PanelResult(p).Rs            = Rs;
    PanelResult(p).Rsh           = Rsh;
    PanelResult(p).nNsVth        = nNsVth;
    PanelResult(p).DiodeResult   = DiodeResult;
    PanelResult(p).PCE           = PCE;
    PanelResult(p).PCE_percent   = PCE_percent;
end

end
