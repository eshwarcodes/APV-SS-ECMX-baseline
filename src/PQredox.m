function [qL,NPQ,Fvp,Fmp,Fop,PHIpsII] = PQredox(J,PAR) %Photosynthesis Research (2019) 141:83–97 
                               % https://doi.org/10.1007/s11120-019-00632-x

% INPUT
% J rate of whole chain electron transport [umol electrons/ s m^2 ]
% PAR absorbed light [umol photons/ s m^2 ]
% OUTPUT
% qL fluorescence parameter indicating redox state of plastoquinone A

% parameters
fPSII = 1.0; % Proportion of absorbed light partitioned to PSII
NPQmax = 2.24;  % Asymptote value sigmoidal fit of non-photochemical quenching (NPQ) response to light intensity
NPQ0 = 0.15; % Basal NPQ value
KNPQ = 1042; % Light intensity at half amplitude of NPQ
nNPQ = 2.52; % Apparent Hill coefficient for NPQ response to light
Fm = 5.0; % dark-adapted maximal fluorescence
Fo = 1.0; % dark-adapted minimal fluorescence
Fv = Fm - Fo; % variable fluorescence (Fv/Fm should be around 0.6-0.8);
m = 2.34E-4; % Slope parameter to estimate effect of reaction centre inactivation on minimal fluorescence (Fp0)
n = 0.038; % Intercept parameter to estimate effect of reaction centre inactivation on Fp0

PHIpsII = J/(PAR*fPSII); % operating efficiency of photosystem II

% NPQ - non-photochemical quenching
if PAR > 0
    NPQ = (NPQmax - NPQ0)/((KNPQ/PAR)^nNPQ + 1) + NPQ0;
else
    NPQ = 0;
end

Fmp = Fm/(NPQ + 1);
Fp = Fmp*(1 - PHIpsII);
FpoNPQ = Fo/(Fv/Fm+Fo/Fmp);
RLDFNPQ = 1-FpoNPQ/Fmp;
RLDF = (1-(m*(0.5*PAR*Fp/Fmp)+n))*RLDFNPQ;
Fvp = RLDF*Fmp;
Fop = Fmp - Fvp;

qL = (Fmp - Fp)*Fop/Fvp/Fp;
% figure(1)
% hold on
% plot(J,Fvp/Fmp,'ko')
return


