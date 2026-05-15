function [dataStruct, hourNumbersSorted] = ...
    loadMatFilesFromDir_range(folder, filePattern, regexPattern, ...
                              startHour, endHour)

files = dir(fullfile(folder, filePattern));

hourNumbers = [];
validFiles  = {};

for i = 1:length(files)

    tokens = regexp(lower(files(i).name), regexPattern, 'tokens');
    if isempty(tokens), continue; end

    h = str2double(tokens{1}{1});

    % Wrap-aware condition
    if startHour <= endHour
        inRange = (h >= startHour && h <= endHour);
    else
        inRange = (h >= startHour || h <= endHour);
    end

    if inRange
        hourNumbers(end+1) = h; %#ok<AGROW>
        validFiles{end+1}  = files(i).name; %#ok<AGROW>
    end
end

if isempty(validFiles)
    warning('No files found in requested range.');
    dataStruct = struct();
    hourNumbersSorted = [];
    return
end

% Sort by hour
[hourNumbersSorted, sortIdx] = sort(hourNumbers);
validFiles = validFiles(sortIdx);

% ------------------------------
% LOAD FILES INTO CELL STRUCTURE
% ------------------------------

dataStruct = struct();

for i = 1:length(validFiles)

    tmp = load(fullfile(folder, validFiles{i}));
    fields = fieldnames(tmp);

    for f = 1:length(fields)

        fieldName = fields{f};

        % Initialize cell array if first time
        if i == 1
            dataStruct.(fieldName) = cell(length(validFiles),1);
        end

        dataStruct.(fieldName){i} = tmp.(fieldName);

    end
end

end