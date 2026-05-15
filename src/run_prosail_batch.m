function [sun_tables, shade_tables] = run_prosail_batch(unique_LAI, zenith1, az_s_matrix1)
    % ProSAIL canopy radiative transfer via Python bridge.
    %
    % Requires:
    %   1. The upstream pyPro4SAIL package on the Python path
    %      (https://github.com/hectornieto/pypro4sail).
    %   2. The thin wrapper prosail_run.py shipped in python/ at the repo
    %      root, which exposes run_prosail_single(LAI, skyl, zen, azim).
    %
    % Resolution order for the wrapper directory:
    %   a) PROSAIL_PY_DIR environment variable (full path to the folder
    %      containing prosail_run.py), if set.
    %   b) <repo>/python (relative to this .m file).
    env_dir = getenv('PROSAIL_PY_DIR');
    if ~isempty(env_dir)
        python_script_path = env_dir;
    else
        here = fileparts(mfilename('fullpath'));
        python_script_path = fullfile(fileparts(here), 'python');
    end
    python_script_path = strrep(python_script_path, '\', '/');

    pyrun("import sys")
    pyrun("sys.path.append('" + python_script_path + "')")
    pyrun("from prosail_run import run_prosail_single")
    % Build angle bins only for simulated hours

    % Initialize result struct
    result = struct;

    for i = 1:length(unique_LAI)
        LAI_val = unique_LAI(i,1);

        % Use the same zenith and azimuth for all entries
        sun_struct = call_prosail(LAI_val, 0.1, zenith1, az_s_matrix1);
        shade_struct = call_prosail(LAI_val, 0.9, zenith1, az_s_matrix1);

        result(i).LAI = LAI_val;
        result(i).sun = sun_struct;
        result(i).shade = shade_struct;
    end

    % Preallocate cell arrays for storing tables
    sun_tables = cell(length(result), 1);
    shade_tables = cell(length(result), 1);

    for i = 1:length(result)
        sun = result(i).sun;
        shade = result(i).shade;

        sun_tables{i} = table( ...
            double(py.array.array('d', py.numpy.nditer(sun{'Wavelength_nm'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Canopy_Reflectance'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Canopy_Transmittance'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Canopy_Absorptance'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Leaf_Reflectance'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Leaf_Transmittance'})))', ...
            double(py.array.array('d', py.numpy.nditer(sun{'Leaf_Absorptance'})))', ...
            'VariableNames', {'Wavelength_nm', 'CanopyReflectance', 'CanopyTransmittance', ...
            'CanopyAbsorptance', 'LeafReflectance', 'LeafTransmittance', 'LeafAbsorptance'});

        shade_tables{i} = table( ...
            double(py.array.array('d', py.numpy.nditer(shade{'Wavelength_nm'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Canopy_Reflectance'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Canopy_Transmittance'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Canopy_Absorptance'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Leaf_Reflectance'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Leaf_Transmittance'})))', ...
            double(py.array.array('d', py.numpy.nditer(shade{'Leaf_Absorptance'})))', ...
            'VariableNames', {'Wavelength_nm', 'CanopyReflectance', 'CanopyTransmittance', ...
            'CanopyAbsorptance', 'LeafReflectance', 'LeafTransmittance', 'LeafAbsorptance'});
    end
end
