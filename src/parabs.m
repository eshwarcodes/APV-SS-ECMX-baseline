function [Q,K2Q,PHI]=parabs(bs,qy,ab) % integrate irradiance spectrum to compute absorbed PAR and leaf available energy

% constant values
h = 6.6261e-34; % Planck constant
c = 3.0e+8; % light speed
Na = 6.02e+23; % Avogadro number

% interpolate bs and a
bsint = fit(bs(:,1),bs(:,2),'linearinterp'); % incident light spectrum [umol photons m-2 s-1 or W m-2]
abint = fit(ab(:,1),ab(:,2),'linearinterp'); % absorptance spectrum [-]
qyint = fit(qy(:,1),qy(:,2),'linearinterp'); % absorptance spectrum [-]

% evaluate functions over entire spectrum (can be expanded toward UV and
% IR)
L = [350:750]';
BS = feval(bsint,L);
AB = feval(abint,L);
QY = feval(qyint,L);

% numerical integration for computation of Q and PHI
ANSGS = 2;
switch ANSGS
    case 1     % from PFD to IRRADIANCE
        F1 = BS.*AB;
        F2 = F1.*QY;
        F3 = 1.0e-6.*F1.*h.*c./(1.0e-9*L).*Na;
        Q = trapz(L,F1); % absorbed PPDF [umol photons m-2 s-1]
        K2Q = trapz(L,F2); % [umol electrons m-2 s-1]
        PHI = trapz(L,F3); % Leaf Available Energy [W m-2]
    case 2     % from IRRADIANCE to PFD
        F1 = BS.*AB;
        F2 = (1.0e-9*L).*F1./(1.0e-6.*h.*c.*Na);
        F3 = F2.*QY;
        PHI = trapz(L,F1); % Leaf Available Energy [W m-2]
        Q = trapz(L,F2); % absorbed PPDF [umol photons m-2 s-1]
        K2Q = trapz(L,F3); % [umol electrons m-2 s-1]
end

return