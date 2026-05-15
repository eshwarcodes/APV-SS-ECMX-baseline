function summary = extract_summary_single(md, Tao_daily, planting_density)
% EXTRACT_SUMMARY_SINGLE  Scalar summary from a single month_data struct.
%
% md                — struct returned by simulateCoupledCropThermal4
% Tao_daily         — [numDays×1] daily mean air temperature (for dTcan)
% planting_density  — plants/m² (lettuce); converts yield from g/plant to kg/m²

    if nargin < 3, planting_density = []; end

    % --- 1) YIELD (last non-empty day) ---
    valid_days = find(~cellfun(@isempty, md.HeadWeight));
    if isempty(valid_days)
        summary = struct( ...
            'meanYield',NaN, 'stdYield',NaN, 'minYield',NaN, ...
            'maxYield',NaN, 'CVYield',NaN, ...
            'harvestMean',NaN, 'harvestStd',NaN, 'harvestIQR',NaN, ...
            'WUE',NaN, 'dTcan',NaN, ...
            'tempMean',NaN, 'tempStd',NaN, 'tempCV',NaN);
        return
    end
    final_grid = md.HeadWeight{valid_days(end)};           % g/plant
    if ~isempty(planting_density)
        final_grid = final_grid * planting_density / 1000; % kg/m2
    end

    summary.meanYield = mean(final_grid(:), 'omitnan');
    summary.stdYield  = std(final_grid(:),  'omitnan');
    summary.minYield  = min(final_grid(:));
    summary.maxYield  = max(final_grid(:));
    if summary.meanYield > 0
        summary.CVYield = summary.stdYield / summary.meanYield;
    else
        summary.CVYield = NaN;
    end

    % --- 2) HARVEST TIMING ---
    hdg = md.harvest_day_grid;
    summary.harvestMean = mean(hdg(:), 'omitnan');
    summary.harvestStd  = std(hdg(:),  'omitnan');
    q = quantile(hdg(:), [0.25 0.75]);
    summary.harvestIQR  = q(2) - q(1);

    % --- 3) WUE ---
    totalYield = sum(final_grid(:), 'omitnan');
    totalET = 0;
    for d = 1:numel(md.ET_mm_day)
        ETgrid = md.ET_mm_day{d};
        if ~isempty(ETgrid)
            totalET = totalET + sum(ETgrid(:), 'omitnan');
        end
    end
    if totalET > 0
        summary.WUE = totalYield / totalET;
    else
        summary.WUE = NaN;
    end

    % --- 4) THERMAL STABILITY ---
    nDays = numel(md.crop_temperature);
    nElem = 0;
    for d = 1:nDays
        Tgrid = md.crop_temperature{d};
        if ~isempty(Tgrid), nElem = nElem + numel(Tgrid); end
    end

    if nElem > 0
        temp_vals = zeros(nElem, 1);
        day_means = zeros(nDays, 1);
        nValid = 0;
        idx = 0;
        for d = 1:nDays
            Tgrid = md.crop_temperature{d};
            if ~isempty(Tgrid)
                n = numel(Tgrid);
                temp_vals(idx+1:idx+n) = Tgrid(:);
                idx = idx + n;
                nValid = nValid + 1;
                day_means(nValid) = mean(Tgrid(:), 'omitnan');
            end
        end
        temp_vals = temp_vals(1:idx);
        day_means = day_means(1:nValid);
        summary.tempMean = mean(temp_vals, 'omitnan');
        summary.tempStd  = std(temp_vals,  'omitnan');
        summary.tempCV   = summary.tempStd / summary.tempMean;

        % dTcan: mean daily canopy temp minus mean daily air temp
        if nargin >= 2 && ~isempty(Tao_daily) && nValid > 0
            summary.dTcan = mean(day_means) - mean(Tao_daily(1:nValid));
        else
            summary.dTcan = NaN;
        end
    else
        summary.tempMean = NaN;
        summary.tempStd  = NaN;
        summary.tempCV   = NaN;
        summary.dTcan    = NaN;
    end

end
