function [PPFD, IPAR, PHI] = parabs_batch(bs_cell, qy, ab)
    % Inputs:
    %   bs_cell - (numY x numX) cell array, each: [λ, I(λ)] (516x2)
    %   qy      - quantum yield spectrum [λ, QY]
    %   ab      - absorption spectrum [λ, Abs]
    % Outputs:
    %   PPFD    - absorbed PPFD [µmol photons/m²/s]
    %   IPAR    - electron flux [µmol e⁻/m²/s]
    %   PHI     - absorbed energy [W/m²]

    [numY, numX] = size(bs_cell);
    PPFD = zeros(numY, numX);
    IPAR = zeros(numY, numX);
    PHI  = zeros(numY, numX);

    % Common wavelength grid [nm]
    L = (350:750)';
    L_m = L * 1e-9;

    % Interpolated absorption and QY
    AB = interp1(ab(:,1), ab(:,2), L, 'linear', 0);
    QY = interp1(qy(:,1), qy(:,2), L, 'linear', 0);

    % Physical constants
    h  = 6.6261e-34;       % J·s
    c  = 3e8;              % m/s
    NA = 6.022e23;         % mol⁻¹

    % Energy per photon [J/photon]
    E_photon = h * c ./ L_m;

    for y = 1:numY
        for x = 1:numX
            bs = bs_cell{y,x};
            if isempty(bs)
                continue;
            end

            % Interpolate spectral irradiance
            I_interp = interp1(bs(:,1), bs(:,2), L, 'linear', 0);  % W/m²/nm

            % Absorbed power spectrum
            A = I_interp .* AB;                 % W/m²/nm

            % Convert to µmol photons/m²/s/nm
            N_photons = (A ./ E_photon) * 1e6 / NA;    % µmol photons/m²/s/nm
            N_electrons = N_photons .* QY;

            % Integrate over wavelength
            PHI(y,x)  = trapz(L, A);             % Absorbed energy [W/m²]
            PPFD(y,x) = trapz(L, N_photons);     % PPFD [µmol photons/m²/s]
            IPAR(y,x) = trapz(L, N_electrons);   % Electron flux [µmol e⁻/m²/s]
        end
    end
end
