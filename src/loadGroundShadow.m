function S = loadGroundShadow(filename, trackingType)
    % loadGroundShadow - Load ground shadow data into a struct.
    % Safe for use inside parfor.
    %
    % Inputs:
    %   filename     - Path to .mat file
    %   trackingType - Tracker type ('hm-panel-system' or other)
    %
    % Output:
    %   S - Struct containing loaded shadow data

    arguments
        filename (1,:) char
        trackingType (1,:) char
    end

    % Base variables common to all trackers
    baseVars = {'panelCenters', 'panelCorners', 'shadowPolygons'};

    if strcmpi(trackingType, 'hm-panel-system')
        extraVars = {'hm1Center', 'hm2Center', 'hmCorners', 'hmShadows', 'entityShadows'};
        varsToLoad = [baseVars, extraVars];
    else
        varsToLoad = baseVars;
    end

    S = load(filename, varsToLoad{:});  % Return as struct
end
