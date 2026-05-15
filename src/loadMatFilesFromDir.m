function [outputStruct,hourNumbersSorted, sortIdx] = loadMatFilesFromDir(dataFolder, filePattern, regexPattern)
% loadMatFilesFromDir  Load and sort .mat files from a folder and aggregate variables.
%
%   outputStruct = loadMatFilesFromDir(dataFolder, filePattern, regexPattern)
%
%   INPUTS:
%       dataFolder  - String specifying the folder where the .mat files reside.
%       filePattern - (Optional) File pattern to match the .mat files.
%                     For example: 'irradianceGround_hour_*.mat'. Default is '*.mat'.
%       regexPattern- (Optional) A regular expression that extracts a numerical
%                     token from the file names used for sorting.
%                     For example: 'irradianceground_hour_(\d+)\.mat'. If omitted or empty,
%                     no sorting is performed.
%
%   OUTPUT:
%       outputStruct - A structure where each field corresponds to a variable loaded from the .mat files.
%                      Each field is a cell array (named as <varname>_all) with one entry per file.
%
%   Example:
%       % Load files and sort by the hour number in the file names:
%       outData = loadMatFilesFromDir(pwd, 'irradianceGround_hour_*.mat', 'irradianceground_hour_(\d+)\.mat');
%
%   Note:
%       This function assumes that every .mat file contains the same set of variables.
%       If a file is missing one of the expected variables, a warning is displayed and an empty
%       value is stored for that file.
%

    % Validate and set defaults for input arguments
    if nargin < 1
        error('Data folder must be provided.');
    end
    if nargin < 2 || isempty(filePattern)
        filePattern = '*.mat';
    end

    % Get list of .mat files in the folder matching filePattern
    files = dir(fullfile(dataFolder, filePattern));
    if isempty(files)
        error('No .mat files found in the folder: %s', dataFolder);
    end

    % If a regex pattern is provided, extract numeric tokens for sorting the files.
    if nargin >= 3 && ~isempty(regexPattern)
        numFiles = length(files);
        sortValues = zeros(numFiles, 1);
        for k = 1:numFiles
            nameLower = lower(files(k).name);  % case-insensitive matching
            tokens = regexp(nameLower, regexPattern, 'tokens');
            if ~isempty(tokens)
                sortValues(k) = str2double(tokens{1}{1});
            else
                error('Filename "%s" does not match expected format defined by the regex.', files(k).name);
            end
        end
        [hourNumbersSorted, sortIdx] = sort(sortValues);
        sortedFiles = files(sortIdx);
    else
        % If no regex is provided, preserve the original order given by dir.
        sortedFiles = files;
    end

    % Load first file to determine which variables to aggregate.
    firstFileName = fullfile(dataFolder, sortedFiles(1).name);
    S_first = load(firstFileName);
    fields = fieldnames(S_first);

    % Preallocate output structure.
    % For each variable in the file, create a cell array in output structure.
    outputStruct = struct();
    for i = 1:length(fields)
        fieldName = fields{i};
        outputStruct.([fieldName, '_all']) = cell(length(sortedFiles), 1);
    end

    % Loop over each file, load the variables, and store them in the output structure.
    for k = 1:length(sortedFiles)
        fullFileName = fullfile(dataFolder, sortedFiles(k).name);
        S = load(fullFileName);

        for i = 1:length(fields)
            fieldName = fields{i};
            if isfield(S, fieldName)
                outputStruct.([fieldName, '_all']){k} = S.(fieldName);
            else
                warning('Variable "%s" not found in file: %s', fieldName, sortedFiles(k).name);
                outputStruct.([fieldName, '_all']){k} = [];  % Optional: store empty if variable not present
            end
        end
    end

end
