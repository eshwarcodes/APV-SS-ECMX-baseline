function resultsPhotosynthesis = runPhotosynthesisGrid(Ri_spec, T_air, ri, qy, ab, ...
    Ca, Csl, ra, rb, Pre, RH, CT, T0, Vcmax0, Oa, g1, g0, rjv, theta, alpha,wind_speed,...
    LAI_t)
    % Constants
    gamma = 0.066;            % psychrometric constant (kPa/°C)
    lambda = 2.45e6;          % latent heat of vaporization (J/kg)
    R_n = 100;                % net radiation at crop surface (W/m²) — placeholder

    [numGroundY, numGroundX] = size(Ri_spec);

    % Preallocate result matrices
    An     = zeros(numGroundY, numGroundX);
    CcF    = zeros(numGroundY, numGroundX);
    Rdark  = zeros(numGroundY, numGroundX);
    gsCO2  = zeros(numGroundY, numGroundX);
    J      = zeros(numGroundY, numGroundX);
    rs     = zeros(numGroundY, numGroundX);
    NPQ    = zeros(numGroundY, numGroundX);
    Fvp    = zeros(numGroundY, numGroundX);
    Fmp    = zeros(numGroundY, numGroundX);
    r_s_pm_data = zeros(numGroundY, numGroundX); 
    r_a_pm_data = zeros(numGroundY, numGroundX);
    Q_ET_data = zeros(numGroundY, numGroundX); 
    transpirationVal    = zeros(numGroundY, numGroundX);
    wlngth = ri(:,1);
    T_crop_K = T_air + 273.15;
    T_amb_K = T_air + 273.15;
    % Loop through grid points (parallel over X)
    for y = 1:numGroundY
        parfor x = 1:numGroundX
            spectrum = Ri_spec{y,x};  % [W/m²/nm]
            bs = [wlngth, spectrum];
            [PPFD_val, IPAR_val, ~] = parabs(bs, qy, ab);

            [Cc, An_val, rs_layer, Rdark_val, gs, J_val, npq, fvp, fmp,transpirationLeaf] = photosynthesis(...
                Ca, IPAR_val, PPFD_val, Csl, ra, rb, T_air, Pre, RH, ...
                CT, T0, Vcmax0, Oa, g1, g0, rjv, theta, alpha);
            [r_s_pm, r_a_pm] = ...
                derivePenmanResistances(gs,...
                ra, rb,...
                T_air, Pre);
            Q_ET = ...
            estimateEvapotranspiration(T_crop_K,...
            T_amb_K, RH, wind_speed,...
            LAI_t,...
            gamma,lambda,R_n,r_s, r_a)

            % Store in matrices
            An(y,x)     = An_val;
            CcF(y,x)    = Cc;
            Rdark(y,x)  = Rdark_val;
            gsCO2(y,x)  = gs;
            J(y,x)      = J_val;
            rs(y,x)     = rs_layer;
            NPQ(y,x)    = npq;
            Fvp(y,x)    = fvp;
            Fmp(y,x)    = fmp;
            transpirationVal(y,x) = transpirationLeaf;
            Q_ET_data(y,x) = Q_ET;
            r_s_pm_data(y,x) = r_s_pm; 
            r_a_pm_data(y,x) = r_a_pm;

        end
    end

    % Pack into output struct
    resultsPhotosynthesis = struct( ...
        'An', An, ...
        'CcF', CcF, ...
        'Rdark', Rdark, ...
        'gsCO2', gsCO2, ...
        'J', J, ...
        'rs', rs, ...
        'NPQ', NPQ, ...
        'Fvp', Fvp, ...
        'Fmp', Fmp, ...
        'transpirationVal',transpirationVal, ...
        'Q_ET_data',Q_ET_data,...
        'r_s_pm_data',r_s_pm_data,...
        'r_a_pm_data',r_a_pm_data);
end
