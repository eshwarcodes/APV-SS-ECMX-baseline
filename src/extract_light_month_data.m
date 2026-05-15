function month_light = extract_light_month_data(month_data)
% Accepts either a single struct or a cell array of structs.

% --- Single-struct path (current simulateCoupledCropThermal4 output) ---
if isstruct(month_data) && ~iscell(month_data)
    month_light = struct();
    month_light.HeadWeight       = month_data.HeadWeight;
    month_light.crop_temperature = month_data.crop_temperature;
    month_light.ET_mm_day        = month_data.ET_mm_day;
    month_light.ET_mm_day_Rn     = month_data.ET_mm_day_Rn;
    return
end

% --- Legacy cell-array path ---
nMonths = length(month_data);
month_light = cell(size(month_data));

for m = 1:nMonths

    if isempty(month_data{m})
        month_light{m} = [];
        continue
    end

    tmp = struct();

    tmp.HeadWeight        = month_data{m}.HeadWeight;
    tmp.crop_temperature  = month_data{m}.crop_temperature;
    tmp.ET_mm_day         = month_data{m}.ET_mm_day;
    tmp.ET_mm_day_Rn      = month_data{m}.ET_mm_day_Rn;

    month_light{m} = tmp;

end

end