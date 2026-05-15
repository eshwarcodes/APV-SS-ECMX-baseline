function S = loadViewFactorGroundPV(filename)
    % loadViewFactorGroundPV - Load view factor ground-PV data into a struct.
    % Safe for use in parfor loops.
    %
    % Input:
    %   filename - full path to viewFactor_hour_X.mat file
    %
    % Output:
    %   S - struct containing sparse view factor variables

    arguments
        filename (1,:) char
    end

    S = load(filename, ...
        'viewFactor_FrontPVGround_sparse', ...
        'viewFactor_RearPVGround_sparse', ...
        'viewFactor_FrontHM1Ground_sparse', ...
        'viewFactor_RearHM1Ground_sparse', ...
        'viewFactor_FrontHM2Ground_sparse', ...
        'viewFactor_RearHM2Ground_sparse');
end