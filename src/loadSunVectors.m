function S = loadSunVectors(filename, trackingType)
    % loadSunVectors - Load sun vector data into a struct (safe for parfor).
    %
    % Inputs:
    %   filename     - Full path to .mat file
    %   trackingType - Tracker type ('hm-panel-system' or other)
    %
    % Output:
    %   S - Struct containing all loaded fields

    arguments
        filename (1,:) char
        trackingType (1,:) char
    end

    % Common variables
    baseVars = {'n_panel', 'sunDir', 'v_shadow', 'panelTilt', 'panelAzimuth'};

    if strcmpi(trackingType, 'hm-panel-system')
        extraVars = {'n_hm1', 'n_hm2', 'hm1Tilt', 'hm2Tilt', 'hm1Azimuth', 'hm2Azimuth'};
        varsToLoad = [baseVars, extraVars];
    else
        varsToLoad = baseVars;
    end

    % Load into struct and return
    S = load(filename, varsToLoad{:});
end
