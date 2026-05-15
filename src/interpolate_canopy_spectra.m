function [canopy_refl_interp, canopy_trans_interp] = interpolate_canopy_spectra(LAI, ref_LAIs, sun_interp)
%INTERPOLATE_CANOPY_SPECTRA Interpolate canopy reflectance and transmittance
%   Inputs:
%       LAI         - [M x N] matrix of LAI values
%       ref_LAIs    - [1 x K] or [K x 1] array of reference LAI values
%       sun_interp  - [1 x K] struct array with fields:
%                       - CanopyReflectance: [B x 1] double
%                       - CanopyTransmittance: [B x 1] double
%
%   Outputs:
%       canopy_refl_interp  - [M x N] cell array, each cell is [B x 1] double
%       canopy_trans_interp - [M x N] cell array, each cell is [B x 1] double

    [m, n] = size(LAI);
    num_LAI = length(ref_LAIs);
    num_bands = length(sun_interp(1).CanopyReflectance);

    % Stack spectra for faster interpolation
    refl_stack = zeros(num_LAI, num_bands);
    trans_stack = zeros(num_LAI, num_bands);

    for i = 1:num_LAI
        refl_stack(i, :) = sun_interp(i).CanopyReflectance;
        trans_stack(i, :) = sun_interp(i).CanopyTransmittance;
    end

    % Prepare output cell arrays
    canopy_refl_interp = cell(m, n);
    canopy_trans_interp = cell(m, n);

    % Interpolate spectra for each LAI value
    for row = 1:m
        for col = 1:n
            lai_val = LAI(row, col);

            refl_interp = interp1(ref_LAIs, refl_stack, lai_val, 'linear', 'extrap');
            trans_interp = interp1(ref_LAIs, trans_stack, lai_val, 'linear', 'extrap');
            refl_interp(refl_interp<0)=0;
            trans_interp(trans_interp<0)=0;
            canopy_refl_interp{row, col} = refl_interp(:);  % Ensure column vector
            canopy_trans_interp{row, col} = trans_interp(:);
        end
    end
end
