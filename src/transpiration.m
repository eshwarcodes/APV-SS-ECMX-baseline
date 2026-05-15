function [T]=transpiration(ra,rb,rs,RH,PHI,Ta,Pa)


%%% INPUT
%%% ra =[s/m] air resistance
%%% rb =[s/m] leaf boundary resistance
%%% rs =[s/m] stomatal resistance
%%% Ta = air temperature [°C]
%%% Pa = atmospheric pressure
%%% RH = relative humidity

%%% OUTPUT
%%% T = [m/s] transpiration rate

%%% PARAMETERS
Rd = 287.058; % [J kg-1 K-1] specific gas constant for dry air
Rv = 461.5; % [J kg-1 K-1] specific gas constant for dry air
lw = 2.5e+6; % [J kg-1] latent heat of water vaporization
cp = 1012; % [J kg-1 K-1] specific heat of air
rho = 1.2; % [kg m-3] air density
rhow = 1000; % [kg m-3] water density

% compute specific humidity from RH
es=0.6108*exp(17.27*Ta/(Ta+237.3))*1000; % saturation vapor pressure [Pa]
ea=RH/100*es; %vapor pressure [Pa]
w=ea*Rd/(Rv*(Pa-ea)); %mixing ratio
q=w/(w+1); % specific humidity
ws=es*Rd/(Rv*(Pa-es)); %mixing ratio at saturation
qs=ws/(ws+1); % specific humidity at saturation
D=qs-q; % vapor pressure deficit in terms of specific humidity [-]

gw=Pa*cp/(0.622*lw); % [Pa K-1] psychrometric constant
S=4098*es/(Ta+237.3)^2; % slope of saturation vapor pressure curve [Pa K-1]
gs=1/rs;
gba=1/(ra+rb);

T=(lw*gw*gba*rho*D+S*PHI)*gs/(rhow*lw*(gw*(gba+gs)+gs*S));
return



