function [harvest_day_grid, HeadWeight_grid, Wg, Ws, GDD, crop_data] = simulateLettuceGrowthSingleDay( ...
    ri, sp, T_air, PAR, ...
    target_head_weight, planting_density, dry_matter_fraction, usable_fraction, ...
    Ca, Csl, Oa, ra, rb, Pre, RH, ...
    CT, T0, Vcmax0, g1, g0, rjv, theta, alpha, ...
    Ws, Wg, GDD, ...
    c_shoot, c_root, f_root, cp, ...
    c_gr_max, Tbase, c_lar)

% Setup
qy = [sp(:,1), sp(:,2)];
ab = [sp(:,1), sp(:,3)];
ri1 = ri(:,1);
ri9 = ri(:,9);

[numGroundY, numGroundX] = size(PAR);
N = numGroundY * numGroundX;

T_air_flat = reshape(T_air, [N, 1]);
PAR_flat = reshape(PAR, [N, 1]);
Wg_flat = reshape(Wg, [N, 1]);
Ws_flat = reshape(Ws, [N, 1]);
GDD_flat = reshape(GDD, [N, 1]);

bs_all = arrayfun(@(par_val) [ri1, (par_val * 1e6 / 86400) * ri9], PAR_flat, 'UniformOutput', false);
[PPFD_vals, IPAR_vals, ~] = parabs_batch(bs_all, qy, ab);
[CcF, An_vec, ~, Rdark_vec, ~, ~, ~, ~, ~] = photosynthesis_vec(...
    Ca * ones(N,1), IPAR_vals, PPFD_vals, Csl * ones(N,1), ra * ones(N,1), rb * ones(N,1), T_air_flat, ...
    Pre * ones(N,1), RH * ones(N,1), CT, T0, Vcmax0, Oa * ones(N,1), g1, g0, rjv, theta, alpha);

An_vec(PPFD_vals < 1e-3) = -Rdark_vec(PPFD_vals < 1e-3);
f_light = 1 - exp(-0.7 * c_lar * Wg_flat);
Pg_vec = An_vec .* 12 ./ 1e6 .* 86400 .* f_light;

Rm = (c_shoot * (1 - f_root) + c_root * f_root) .* (Wg_flat + Ws_flat) .* 2.0.^((T_air_flat - 25)/10);
gr_temp = 1.6.^((T_air_flat - 20)/10);
frac = Ws_flat ./ (Wg_flat + Ws_flat);
r_gr = c_gr_max * frac .* gr_temp * 86400;
dWg = min(r_gr .* Ws_flat, Ws_flat);

Wg_flat = Wg_flat + dWg;
Ws_flat = Ws_flat + (Pg_vec - cp * dWg - Rm);
GDD_flat = GDD_flat + max(0, T_air_flat - Tbase);

Biomass_fresh = Wg_flat ./ dry_matter_fraction * usable_fraction;
head_weight = Biomass_fresh ./ planting_density;
harvest_day = double(head_weight >= target_head_weight);
harvest_day(harvest_day == 0) = NaN;

% Reshape back to grid
Pg_grid = reshape(Pg_vec, numGroundY, numGroundX);
HeadWeight_grid = reshape(head_weight, numGroundY, numGroundX);
Wg = reshape(Wg_flat, numGroundY, numGroundX);
Ws = reshape(Ws_flat, numGroundY, numGroundX);
GDD = reshape(GDD_flat, numGroundY, numGroundX);
harvest_day_grid = reshape(harvest_day, numGroundY, numGroundX);

crop_data = struct(...
    'Pg', Pg_grid, 'Wg', Wg, 'Ws', Ws, ...
    'HeadWeight', HeadWeight_grid, 'GDD', GDD, ...
    'harvested', harvest_day_grid);
end
