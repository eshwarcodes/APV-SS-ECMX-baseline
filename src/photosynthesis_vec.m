% Vectorized Photosynthesis Function with Case Options
function [CcF, An, rs, Rdark, gsCO2, J, NPQ, Fvp, Fmp, E] =...
    photosynthesis_vec(Cc, IPAR, PPFD, Csl, ra, rb, Ts, Pre, RH, CT, T0, Vmax, Oa, g1, go, rjv, theta, alpha)

% Vectorization setup
vars = {Cc, IPAR, PPFD, Csl, ra, rb, Ts, Pre, RH};
for i = 1:length(vars)
    if isrow(vars{i}), vars{i} = vars{i}'; end
end
[Cc, IPAR, PPFD, Csl, ra, rb, Ts, Pre, RH] = deal(vars{:});

n = length(Ts);
Tf = 273.15;
Pre0 = 101325;
R = 0.008314;

% Configurable model switches
ANS_TEMP = 1; % 1 = Kattge & Knorr, else Bernacchi
ANSG = 1;     % 0 = Cox, 1 = Leuning, 2 = Bonan
ANSMOD = 2;   % 1 = Fatichi, 2 = Daly
ANSTPU = 1;   % 1 = TPU limited, 2 = no TPU limit
ANSGS = 2;    % 1 = empirical, 2 = qL-based

% Vapor pressure and deficit
es = 0.6108 .* exp(17.27 .* Ts ./ (Ts + 237.3));
Ds = es .* (1 - RH ./ 100);

% Unit conversions
ra = ra .* (0.0224 .* (Ts + Tf) .* Pre0 ./ (Tf .* Pre)) .* 1e-6;
rb = rb .* (0.0224 .* (Ts + Tf) .* Pre0 ./ (Tf .* Pre)) .* 1e-6;
Cc = Cc .* 1e-6 .* Pre;
Csl = Csl .* 1e-6 .* Pre;
Oa = Oa .* 1e-6 .* Pre;
rmes = 1 ./ (1e6 * 0.60);
go = go * 1e6;

% Temperature dependencies
Ts_k = Ts + Tf;
Tref = T0 + Tf;

switch ANS_TEMP
    case 1
        Hd = 200; Ha = 72; DS = 0.650;
    otherwise
        Hd = 149; Ha = 65.33; DS = 0.485;
end
kT_vm = exp(Ha .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k)) .* ...
        (1 + exp((Tref .* DS - Hd) ./ (Tref .* R))) ./ ...
        (1 + exp((Ts_k .* DS - Hd) ./ (Ts_k .* R)));
Vm = Vmax .* kT_vm;

switch ANS_TEMP
    case 1
        Hd = 200; Ha = 50; DS = 0.650;
    otherwise
        Hd = 152; Ha = 43.5; DS = 0.495;
end
kT_jm = exp(Ha .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k)) .* ...
        (1 + exp((Tref .* DS - Hd) ./ (Tref .* R))) ./ ...
        (1 + exp((Ts_k .* DS - Hd) ./ (Ts_k .* R)));
Jmax = Vmax .* rjv;
Jm = Jmax .* kT_jm;

% TPU
TPU25 = 0.1182 .* Vmax;
Ha = 53.1; DS = 0.490; Hd = 150.65;
kT_tpu = exp(Ha .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k)) .* ...
         (1 + exp((Tref .* DS - Hd) ./ (Tref .* R))) ./ ...
         (1 + exp((Ts_k .* DS - Hd) ./ (Ts_k .* R)));
TPU = TPU25 .* kT_tpu;

% CO2 compensation point
switch ANSG
    case 0
        fT = 0.57 .^ (0.1 .* (Ts - 25));
        GAM = (Oa ./ (2 .* 2600 .* fT)) .* (CT == 3);
    case 1
        G0 = 34.6; G1 = 0.0451; G2 = 0.000347;
        GAM = G0 .* (1 + G1 .* (Ts_k - Tref) + G2 .* (Ts_k - Tref).^2);
        GAM = GAM .* 1e-6 .* Pre;
    case 2
        Ha = 37.83;
        kT = exp(Ha .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k));
        GAM25 = 42.75 .* 1e-6 .* Pre;
        GAM = GAM25 .* kT;
end

% Michaelis-Menten constants
Ha_kc = 59.43; Ha_ko = 36.00;
Kc25 = 302 .* 1e-6 .* Pre;
Ko25 = 256 .* 1e-3 .* Pre;
kT_kc = exp(Ha_kc .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k));
kT_ko = exp(Ha_ko .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k));
Kc = Kc25 .* kT_kc;
Ko = Ko25 .* kT_ko;

% Light response model
switch ANSMOD
    case 1
        Q = 0.081 .* IPAR;
        d2 = -(Q + Jm ./ 4); d3 = Q .* Jm ./ 4;
    case 2
        Q = PPFD;
        d2 = -(alpha .* Q + Jm); d3 = alpha .* Q .* Jm;
end
d1 = theta;
disc = sqrt(d2.^2 - 4 .* d1 .* d3);
J = min((-d2 + disc) ./ (2 .* d1), (-d2 - disc) ./ (2 .* d1));

% Rubisco, light, export limitations
JC = Vm .* (Cc - GAM) ./ (Cc + Kc .* (1 + Oa ./ Ko));
if ANSMOD == 2
    JL = (J ./ 4) .* (Cc - GAM) ./ (Cc + 2 .* GAM);
else
    JL = J .* (Cc - GAM) ./ (Cc + 2 .* GAM);
end
JE = 3 .* TPU;

% First polynomial
b1 = 0.9; b2 = -(JC + JL); b3 = JC .* JL;
disc2 = sqrt(b2.^2 - 4 .* b1 .* b3);
JP = min((-b2 + disc2) ./ (2 .* b1), (-b2 - disc2) ./ (2 .* b1));

% Second polynomial
switch ANSTPU
    case 1
        c1 = 0.9; c2 = -(JP + JE); c3 = JP .* JE;
        disc3 = sqrt(c2.^2 - 4 .* c1 .* c3);
        A = min((-c2 + disc3) ./ (2 .* c1), (-c2 - disc3) ./ (2 .* c1));
    case 2
        A = JP;
end

A = A * 1;
% Dark Respiration
if CT == 3
    Ha = 46.39; DS = 0.490; Hd = 150.65;
    Rdark25 = 0.015 .* Vmax;
    kT = exp(Ha .* (Ts_k - Tref) ./ (Tref .* R .* Ts_k)) .* ...
         (1 + exp((Tref .* DS - Hd) ./ (Tref .* R))) ./ ...
         (1 + exp((Ts_k .* DS - Hd) ./ (Ts_k .* R)));
    Rdark = Rdark25 .* kT;
elseif CT == 4
    fT = 2.0 .^ (0.1 .* (Ts - 25));
    fT3 = 1 ./ (1 + exp(1.3 .* (Ts - 55)));
    Rdark25 = 0.025 .* Vmax;
    Rdark = Rdark25 .* fT .* fT3;
end

An = A - Rdark;

switch ANSGS
    case 1
        gsCO2 = go + (1 + g1 ./ sqrt(Ds)) .* An .* Pre ./ Csl;
        NPQ = NaN(n,1); Fvp = NaN(n,1); Fmp = NaN(n,1);
    case 2
        [qL, NPQ, Fvp, Fmp, ~, ~] = arrayfun(@PQredox, J, IPAR);
        gsCO2 = go + (1 + g1 ./ sqrt(Ds)) .* (1 - qL) .* Pre ./ Csl;
end
gsCO2 = max(gsCO2, go);
rsCO2 = 1 ./ gsCO2;
CcF = Csl - An .* Pre .* (rsCO2 + rmes + 1.37 .* rb + ra);
CcF = max(CcF ./ (Pre .* 1e-6), 0);
An = (Csl - CcF .* 1e-6 .* Pre) ./ (Pre .* (rsCO2 + rmes + 1.37 .* rb + ra));

rsH2O = (rsCO2 ./ 1.64) .* 1e6;
rs = rsH2O .* (Tf .* Pre) ./ (0.0224 .* (Ts + Tf) .* Pre0);
gsH2O = (gsCO2 ./ 1e6) .* 1.6;
ea = es .* RH ./ 100;
VPD = max(es - ea, 0);
Pre_kPa = Pre ./ 1000;
E = (VPD .* gsH2O) ./ Pre_kPa;

end
