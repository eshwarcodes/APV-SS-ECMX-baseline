function [Wnow, WFnow, LAInow, Nnow, Ns, Nl, Nf,AgeL, AgeF, Rc_prev,Cstore] = ...
    simulateTomatoGrowthSingleDay_vectorized( ...
        T_canopy, DLI_inside, ...
        Wprev, WFprev, LAIprev, Nprev, ...
        Ns, Nl, Nf,AgeL, AgeF, Rc_prev, day, P, Thours,Smax,eta_store,...
                    eta_withd,Cstore)

% Vectorized tomato day-step over the whole grid.
% Inputs (Y×X unless noted): T_canopy, DLI_inside, Wprev, WFprev, LAIprev, Nprev
% Cohorts: Ns/Nl/Nf are Y×X×nL
% Params P fields used: PD, Nm, LAImax, E, D, LFmax, K, m, Qe, Q10, RL, RF,
%   Tcrit, NFF, CO2, sunrise, sunset, day1, Thours, SLAmin, SLAmax, betaT, betaC

% ---- Sizes ----
[Y, X] = size(Wprev);
nL = size(Ns, 3);

% ---- Constants (scalars) ----
FPN    = P.FPN;
FRPET  = P.FRPET;
FRSTM  = P.FRSTM;
TPL    = P.TPL;
SLAmin = P.SLAmin;  SLAmax = P.SLAmax;  betaT = P.betaT;  betaC = P.betaC;
a1=P.a1; b1=P.b1; c1=P.c1; a2=P.a2; b2=P.b2; c2=P.c2;
tiny = 1e-12;
Tref = P.Tref;
% ---- Unpack params ----
pd     = P.PD;      Nm     = P.Nm;      LAImax = P.LAImax;
E      = P.E;       LFmax  = P.LFmax;
Kpar   = P.K;       m      = P.m;       Qe     = P.Qe;
Q10    = P.Q10;     RL     = P.RL;      RF     = P.RF;
Tcrit  = P.Tcrit;   NFF    = P.NFF;     CO2ppm = P.CO2;
sunrise = P.sunrise;  sunset = P.sunset;  
%sunrise = 5.43;
%sunset = 17.74;


% Map clock hours (possibly fractional) to 1..24 indices on the 24h slice
s_idx = max(1, min(24, floor(sunrise) + 1));   % e.g., 5.8 -> 6
e_idx = max(1, min(24, ceil(sunset)));         % e.g., 18.3 -> 18 (cap at 24)
if e_idx < s_idx, e_idx = s_idx; end           % guard

% Daylength & PPFD
daylen_h = max(sunset - sunrise, 1/60);        % hours (≥ 1 min)
%DLI_inside = 61.13*ones(105,105);
PPFD = DLI_inside * 1e6 / (daylen_h * 3600);   % μmol photons m^-2 s^-1 (Y×X)

% Temps
T1   = T_canopy;                 % canopy temp field (°C)
%T1=14.51*ones(105,105);
Tday = mean(Thours(s_idx:e_idx));% scalar daytime mean for PGREDT
%Tday = 11;% scalar daytime mean for PGREDT

% ---- SLA (m^2 g^-1), bounded ----
SLA = (SLAmin + (SLAmax - SLAmin) .* exp(-0.471 .* DLI_inside)) ./ ...
      ( (1 + betaT.*(24 - T1)) .* (1 + betaC.*(CO2ppm - 350)) );
SLA = max(SLAmin, min(SLA, SLAmax));

% ---- CO2 & temperature factors ----
FCO2_1 = 1 + 0.0003*(CO2ppm - 350);  % scalar
fnT2   = 1 + 0.0281*(T1 - 28);       % scalar
%a1_1=0.025;
%b1_1 = 0.35;
%c1_1=4;
POL3 = a1*b1.*exp(-b1*(AgeL - c1)) .* exp(-exp(-b1*(AgeL - c1)));
% a2_1 = 6;
% b2_1 = 0.1;
% c2_1 = 10;
POF3 = a2*b2.*exp(-b2*(AgeF - c2)) .* exp(-exp(-b2*(AgeF - c2)));

% ---- gTdaytime (temperature reduction factor) ----
% Piecewise response: base ~8 °C, optimum ~25 °C, ceiling ~40 °C
% ---- rF, rL (piecewise) ---- 
rF = zeros(Y,X); 
rL = zeros(Y,X); 
maskCold = T1 < 12; 
maskMid = T1 > 12 & T1 <= Tcrit; 
maskHot = T1 > Tcrit; 
rF(maskCold)=0; 
rL(maskCold)=0; 
rL(maskMid) = (0.0005*T1(maskMid)) - 0.0014; 
rF(maskMid) = (0.0017*T1(maskMid)) - 0.0147; 
rL(maskHot) = max(0, -0.0009*T1(maskHot) + 0.0433); 
rF(maskHot) = max(0, -0.0021*T1(maskHot) + 0.1067); 
rF = min(rF, 0.033); rL = min(rL, 0.013); 
% ---- fnT1 (piecewise) ---- 
fnT1 = zeros(Y,X); 
fnT1(T1<=12) = 0;
if(T1<=12)
    stophere=1;
end
maskA = T1>12 & T1<=28; 
fnT1(maskA) = 1 + 0.0281*(T1(maskA)-28); 
maskB = T1>28; 
fnT1(maskB) = max(0, 1 - 0.0455*(T1(maskB)-28)); 
% ---- PGREDT1 for daytime temperature ---- 
% if Tday < 10 || Tday > 45 
%     PGREDT1 = 0; 
% else 
%     PGREDT1 = -0.0022*Tday^2 + 0.1182*Tday - 0.6401; 
% end

% ---- PGREDT1 based on CANOPY temperature (grid, not ambient) ----
PGREDT1 = zeros(Y,X);
maskPg = (T1 >= 10) & (T1 <= 45);
PGREDT1(maskPg) = -0.0022.*T1(maskPg).^2 + 0.1182.*T1(maskPg) - 0.6401;  % dimensionless
% outside 10..45 °C PGREDT1 stays 0


% ---- fRN1 from nodenumber ----
nodenumber = P.nodeNumber;
fRN        = P.fN;
fRN1 = interp1(nodenumber, fRN, Nprev, 'linear', 'extrap');
fRN1 = max(min(fRN1, max(fRN)), min(fRN));

% ---- Node update ----
dNdt = Nm .* fnT1 * FCO2_1;
Nnow = min(Nprev + dNdt, NFF);

% ---- Flowering gate ----
gateN = Nprev >= NFF;             % Y×X (logical)
gateN3 = repmat(gateN, 1,1,nL);   % expand to Y×X×nL

% ---- POL/POF (expand to 3D for cohort ops) ----
%POLm = POL .* (fnT1 * FCO2_1);     % Y×X
%POFm = POF .* (fnT1 * FCO2_1) .* gateN;
%POL3 = repmat(POL3, 1,1,nL);      % Y×X×nL
%POF3 = repmat(POF3, 1,1,nL);      % Y×X×nL


% Inputs:
%   Ns, Nl, Nf : Y×X×nL double
%   Nprev      : Y×X
%   fnT1, rL, rF : Y×X
%   Nm, FCO2_1, pd, FPN, TPL, NFF, nL : scalars (scalars auto-broadcast)
% Outputs: Ns, Nl, Nf updated in-place

maskFlower = (Nprev >= NFF);           % Y×X logical
% ---- Cohort transfer rates (drop the .*nL multiplier) ----
rLfac = rL .* FCO2_1;    % per day
rFfac = rF .* FCO2_1;    % per day

% ---- Prepare new arrays for today ----
Nl_new   = zeros(Y,X,nL);  AgeL_new = zeros(Y,X,nL);
Nf_new   = zeros(Y,X,nL);  AgeF_new = zeros(Y,X,nL);
Ns_new   = zeros(Y,X,nL);

maskFlower = (Nprev >= NFF);  % Y×X
% ---- Node update & births ----
birthNs = Nm .* fnT1 .* FCO2_1 .* pd;             % stems per m2
birthNl = Nm .* fnT1 .* FCO2_1 .* pd ./ (1+TPL);  % leaves per m2
birthNf = Nm .* fnT1 .* FCO2_1 .* pd .* FPN;      % fruits per m2

for i = 1:nL
    % ---------- STEMS (no age tracked) ----------
    if i==1
        Ns_new(:,:,i) = Ns(:,:,i) + birthNs - rLfac.*Ns(:,:,i);
    elseif i<=nL-1
        Ns_new(:,:,i) = Ns(:,:,i) + rLfac.*Ns(:,:,i-1) - rLfac.*Ns(:,:,i);
    else
        Ns_new(:,:,i) = Ns(:,:,i) + rLfac.*Ns(:,:,i-1);
    end

    % ---------- LEAVES (mass + age) ----------
    if i==1
        birth_m = birthNl;   birth_a = 0;
        inflow_m = 0;        inflow_a = 0;
        stay_m   = (1 - rLfac).*Nl(:,:,i);
        stay_a   = AgeL(:,:,i) + 1;           % age increment happens here
    elseif i<=nL-1
        birth_m = 0;         birth_a = 0;
        inflow_m = rLfac.*Nl(:,:,i-1);
        inflow_a = AgeL(:,:,i-1) + 1;         % incoming cohort already aged +1
        stay_m   = (1 - rLfac).*Nl(:,:,i);
        stay_a   = AgeL(:,:,i) + 1;
    else
        birth_m = 0;         birth_a = 0;
        inflow_m = rLfac.*Nl(:,:,i-1);
        inflow_a = AgeL(:,:,i-1) + 1;
        stay_m   = Nl(:,:,i);                  % last bin no outflow
        stay_a   = AgeL(:,:,i) + 1;
    end

    Nl_new(:,:,i) = birth_m + inflow_m + stay_m;
    denomL = max(Nl_new(:,:,i), 1e-12);
    AgeL_new(:,:,i) = (birth_m.*birth_a + inflow_m.*inflow_a + stay_m.*stay_a) ./ denomL;

    % ---------- FRUITS (mass + age; gate by flowering) ----------
    if i==1
        birthF_m = birthNf .* maskFlower;   birthF_a = 0;
        inflowF_m = 0;                      inflowF_a = 0;
        stayF_m   = (1 - rFfac).*Nf(:,:,i) .* maskFlower;
        stayF_a   = (AgeF(:,:,i) + 1) .* maskFlower;
    elseif i<=nL-1
        birthF_m = 0;                       birthF_a = 0;
        inflowF_m = rFfac.*Nf(:,:,i-1) .* maskFlower;
        inflowF_a = AgeF(:,:,i-1) + 1;
        stayF_m   = (1 - rFfac).*Nf(:,:,i) .* maskFlower;
        stayF_a   = (AgeF(:,:,i) + 1) .* maskFlower;
    else
        birthF_m = 0;                       birthF_a = 0;
        inflowF_m = rFfac.*Nf(:,:,i-1) .* maskFlower;
        inflowF_a = AgeF(:,:,i-1) + 1;
        stayF_m   = (Nf(:,:,i) .* maskFlower);      % last bin no outflow
        stayF_a   = (AgeF(:,:,i) + 1) .* maskFlower;
    end

    Nf_new(:,:,i) = birthF_m + inflowF_m + stayF_m;
    denomF = max(Nf_new(:,:,i), 1e-12);
    AgeF_new(:,:,i) = (birthF_m.*birthF_a + inflowF_m.*inflowF_a + stayF_m.*stayF_a) ./ denomF;
end

% Zero fruits where not flowering
Nf_new(~maskFlower) = 0;  AgeF_new(~maskFlower) = 0;

% Swap in
Ns  = max(Ns_new,0);
Nl  = max(Nl_new,0);  AgeL = max(AgeL_new,0);
Nf  = max(Nf_new,0);  AgeF = max(AgeF_new,0);


% ---- Demands (sum over cohorts) ----
dAlpdt      = Nl .* POL3;                                % Y×X×nL
SLA3        = repmat(max(SLA, tiny), 1,1,nL);                % broadcast SLA
Ldem_matrix = (1+FRPET) .* (dAlpdt ./ SLA3);
Sdem_matrix = FRSTM .* Ldem_matrix .* (Ns ./ max(Nl, tiny));
Fdem_matrix = POF3 .* Nf;

Ldem = sum(Ldem_matrix, 3);      % Y×X
Sdem = sum(Sdem_matrix, 3);
Fdem = sum(Fdem_matrix, 3);
Demand = Ldem + Sdem + Fdem;

% ---- Photosynthesis & respiration (g m^-2 d^-1) ----
Pg = ((sunset-sunrise) * LFmax * PGREDT1 / Kpar) .* ...
     log((((1-m)*LFmax*(sunset-sunrise)) + (Qe*(sunset-sunrise)*Kpar*PPFD)) ./ ...
          (((1-m)*LFmax*(sunset-sunrise)) + (Qe*(sunset-sunrise)*Kpar*PPFD .* exp(-Kpar*LAIprev))));

% Use explicit reference 20°C (adjust if you want 2°C like original)

Rm = (Q10 .^ ((T1 - Tref)/10)) .* (RL .* Wprev + RF .* WFprev);


% ================== STORAGE + ALLOCATION (unit-consistent) ==================
tiny  = 1e-9;
tiny2 = 1e-12;

% ---- Maintenance first (carbon units) ----
balance = Pg - Rm;                                    % carbon (g m^-2 d^-1)
withdraw_rm = min(max(-balance,0), Cstore) .* eta_withd;  % carbon withdrawn for Rm
Cstore = Cstore - withdraw_rm;

Seff   = balance + withdraw_rm;                       % carbon after Rm
GRnet  = E .* max(Seff, 0);                           % growth units available from Pg (+Rm storage)

% ---- Rc diagnostic (nonnegative; just for logging) ----
Demand_eff = max(Demand, tiny);                       % growth units
Rc_prev = GRnet ./ Demand_eff;                        % >= 0

% ---- Fruit mask ----
hasFruit = (Fdem > tiny2);
noFruit  = ~hasFruit;

% ---- Veg allocation with GRnet budget (sink-limited); allow storage top-up per cell if GRnet≈0 ----
Dveg = Ldem + Sdem;                                   % growth units
Rc_LS_noF = min(GRnet ./ max(Dveg, tiny2), 1);        % base veg-satisfaction (no fruit case)
Rc_LS_F   = Rc_LS_noF;                                % same base for fruit case

% Per-cell storage boost used for VEG growth (growth units), tracked for leftover calc
boost_LS = zeros(size(GRnet));                        % growth units added from storage to veg

% --- Cells with NO FRUIT and GRnet ~ 0: try funding veg from storage ---
mask_noF   = (Dveg > tiny2) & (GRnet <= tiny2) & (~hasFruit);
if any(mask_noF(:))
    % available growth from storage today (carbon -> growth via E*eta_withd)
    avail_grow = eta_withd .* E .* Cstore(mask_noF);  % growth units
    need_grow  = Dveg(mask_noF);                      % growth units
    take_grow  = min(avail_grow, need_grow);          % growth units

    % consume storage carbon accordingly
    Cstore(mask_noF) = Cstore(mask_noF) - (take_grow ./ (E*eta_withd));  % carbon units
    Rc_LS_noF(mask_noF) = take_grow ./ max(need_grow, tiny2);            % 0..1
    boost_LS(mask_noF)  = take_grow;                                     % track boost
end

% --- Cells WITH FRUIT and GRnet ~ 0: fund veg from storage so fruit can get any leftover ---
mask_withF_zero = (Dveg > tiny2) & (GRnet <= tiny2) & hasFruit;
if any(mask_withF_zero(:))
    avail_grow = eta_withd .* E .* Cstore(mask_withF_zero);  % growth units
    need_grow  = Dveg(mask_withF_zero);                      % growth units
    take_grow  = min(avail_grow, need_grow);                 % growth units

    Cstore(mask_withF_zero) = Cstore(mask_withF_zero) - (take_grow ./ (E*eta_withd)); % carbon
    Rc_LS_F(mask_withF_zero) = take_grow ./ max(need_grow, tiny2);                    % 0..1
    boost_LS(mask_withF_zero) = take_grow;                                           % growth units
end

% ---- Case 1: no fruit -> veg only, limited by (GRnet or storage-boosted) ----
Lalloc_noF = Ldem .* Rc_LS_noF;
Salloc_noF = Sdem .* Rc_LS_noF;
Falloc_noF = zeros(size(GRnet));
% zero-out where fruit exists
Lalloc_noF(hasFruit) = 0; Salloc_noF(hasFruit) = 0;

% ---- Case 2: fruit present -> satisfy veg sinks, leftover goes to fruit ----
Lalloc_F = Ldem .* Rc_LS_F;
Salloc_F = Sdem .* Rc_LS_F;

% Leftover growth budget for fruit = GRnet + any veg-boost from storage − veg allocations
% (veg-boost was already paid from storage above)
GR_budget_plus = GRnet + boost_LS;                     % total growth attainable today
leftover_F = max(GR_budget_plus - (Lalloc_F + Salloc_F), 0);   % growth units
Falloc_F   = leftover_F;                               % optional: cap with Fdem -> min(leftover_F, Fdem)

% zero-out where no fruit
Lalloc_F(noFruit) = 0; Salloc_F(noFruit) = 0; Falloc_F(noFruit) = 0;

% ---- Combine allocations ----
Lalloc = Lalloc_noF + Lalloc_F;                        % growth units
Salloc = Salloc_noF + Salloc_F;
Falloc = Falloc_noF + Falloc_F;

% ---- Senescence gating on vegetative (if desired) ----
Lalloc = Lalloc .* (1 - fRN1);
Salloc = Salloc .* (1 - fRN1);

% ---- Carbon bookkeeping after allocations ----
used_today   = Lalloc + Salloc + Falloc;               % growth units actually realized
surplus_grow = max(GR_budget_plus - used_today, 0);    % unused growth budget
surplus_car  = surplus_grow ./ max(E, tiny);           % convert to carbon for storage
Cstore       = min(Cstore + eta_store .* surplus_car, Smax);

% ---- State updates ----
LAInow = min(max(LAIprev + Lalloc .* SLA, 0), LAImax);
WFnow  = max(WFprev + Falloc, 0);
Wnow   = max(Wprev  + Lalloc + Salloc, 0);
% ===========================================================================



if(day==80)
    stophere=1;
end

% Numerical guards
LAInow = max(LAInow, 0);
WFnow  = max(WFnow,  0);
Wnow   = max(Wnow,   0);
Nnow   = max(Nnow,   0);
Ns     = max(Ns,     0);
Nl     = max(Nl,     0);
Nf     = max(Nf,     0);

end
