function S = loadSkyViewFactor(filename, trackingType)
    % loadSkyViewFactor - Load SVF variables into a struct (safe for parfor).
    %
    % Inputs:
    %   filename     - Full path to the .mat file
    %   trackingType - 'hm-panel-system' or other
    %
    % Output:
    %   S - Struct containing loaded SVF fields

    arguments
        filename (1,:) char
        trackingType (1,:) char
    end

    baseVars = {'SVF_front', 'SVF_back'};

    if strcmpi(trackingType, 'hm-panel-system')
        extraVars = {'SVF_HM1front', 'SVF_HM1back', 'SVF_HM2front', 'SVF_HM2back'};
        varsToLoad = [baseVars, extraVars];
    else
        varsToLoad = baseVars;
    end

    S = load(filename, varsToLoad{:});
end
