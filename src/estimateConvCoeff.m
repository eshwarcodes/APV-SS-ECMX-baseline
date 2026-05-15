function h = estimateConvCoeff(U, Length, Width,Tao, Pr, numPanels)
  for i=1:numPanels
        L = max(Length, Width); % Characteristic length assumption (m)                 
        k_air = 0.024 + 7e-5 * (Tao(i,1));         % Thermal conductivity [W/m·K]
        nu_air = 1.326e-5 + 9e-8 * (Tao(i,1));     % Kinematic viscosity [m²/s]

        Re = U * L / nu_air;
        if Re < 5e5
            Nu = 0.664 * Re^0.5 * Pr^(1/3);  % laminar
        else
            Nu = 0.037 * Re^0.8 * Pr^(1/3);  % turbulent
        end
        h(i,1) = k_air / L * Nu;
  end
end