function summary = extract_summary(month_data)
% ==========================================================
% EXTRACT_SUMMARY
% Accepts either:
%   - a single struct (sequential day indexing, from simulateCoupledCropThermal4)
%   - a cell array of structs (legacy monthly format)
% ==========================================================

% --- Single-struct path (current simulateCoupledCropThermal4 output) ---
if isstruct(month_data) && ~iscell(month_data)
    summary = extract_summary_single(month_data);
    return
end

% --- Legacy cell-array path ---
nMonths = length(month_data);
summary = struct();

for m = 1:nMonths

    if isempty(month_data{m})
        continue
    end
    
    % ======================================================
    % 1) FINAL YIELD (last available day of month)
    % ======================================================
    
    head_cells = month_data{m}.HeadWeight;
    
    if iscell(head_cells)
        valid_days = find(~cellfun(@isempty, head_cells));
        if isempty(valid_days)
            continue
        end
        final_day = valid_days(end);
        final_grid = head_cells{final_day};
    else
        final_grid = head_cells;
    end
    
    summary(m).meanYield = mean(final_grid(:), 'omitnan');
    summary(m).stdYield  = std(final_grid(:),  'omitnan');
    summary(m).minYield  = min(final_grid(:));
    summary(m).maxYield  = max(final_grid(:));
    
    if summary(m).meanYield > 0
        summary(m).CVYield = summary(m).stdYield / summary(m).meanYield;
    else
        summary(m).CVYield = NaN;
    end
    
    % ======================================================
    % 2) HARVEST TIMING
    % ======================================================
    
    harvest_grid = month_data{m}.harvest_day_grid;
    
    summary(m).harvestMean = mean(harvest_grid(:), 'omitnan');
    summary(m).harvestStd  = std(harvest_grid(:),  'omitnan');
    
    q = quantile(harvest_grid(:), [0.25 0.75]);
    summary(m).harvestIQR = q(2) - q(1);
    
    % ======================================================
    % 3) WATER USE EFFICIENCY
    % ======================================================
    
    totalYield = sum(final_grid(:), 'omitnan');
    totalET = 0;
    
    ET_cells = month_data{m}.ET_mm_day;
    
    for d = 1:length(ET_cells)
        ETgrid = ET_cells{d};
        if isempty(ETgrid), continue; end
        totalET = totalET + sum(ETgrid(:), 'omitnan');
    end
    
    if totalET > 0
        summary(m).WUE = totalYield / totalET;
    else
        summary(m).WUE = NaN;
    end
    
    % ======================================================
    % 4) THERMAL STABILITY
    % ======================================================
    
    temp_cells = month_data{m}.crop_temperature;
    temp_vals = [];
    
    for d = 1:length(temp_cells)
        Tgrid = temp_cells{d};
        if isempty(Tgrid), continue; end
        temp_vals = [temp_vals; Tgrid(:)];
    end
    
    if ~isempty(temp_vals)
        summary(m).tempMean = mean(temp_vals, 'omitnan');
        summary(m).tempStd  = std(temp_vals,  'omitnan');
        summary(m).tempCV   = summary(m).tempStd / summary(m).tempMean;
    else
        summary(m).tempMean = NaN;
        summary(m).tempStd  = NaN;
        summary(m).tempCV   = NaN;
    end
    
end

end