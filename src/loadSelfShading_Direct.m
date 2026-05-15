function S = loadSelfShading_Direct(filename, trackingType)
    % loadSelfShading_Direct - Load self-shading fractions into a struct.
    % Safe for use in parfor loops.

    arguments
        filename (1,:) char
        trackingType (1,:) char
    end

    if strcmpi(trackingType, 'hm-panel-system')
        S = load(filename, 'shadingFractions_direct', ...
                         'shadingFractionsHM1_direct', ...
                         'shadingFractionsHM2_direct');
    else
        S = load(filename, 'shadingFractions_direct');
        S.shadingFractionsHM1_direct = [];  % Pad for consistency
        S.shadingFractionsHM2_direct = [];
    end
end
