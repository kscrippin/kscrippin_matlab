 function [Tmat,Pmat, Ttmat, Ptmat, Mmat, Vmat, Emat] = TFPCA(M0,Alt, Tt4, PIc, PIf, alpha)
%[Tmat,Pmat, Ttmat, Ptmat, Mmat, Vmat, Emat] = TFPCA(M0,Alt, Tt4, PIc, PIf, alpha)
%[degR, psia, degR,psia, unitless, ft/s]=fcn(unitless, ft, degR,unitless, unitless, unitless)

%% Constants
R = 53.35; %(lbf ft)/(lbm degR)
cp = 0.24; % Btu/(lbm ft/s^2)
gamma = 1.4; %unitless
Hpr = 18400; %Btu/lbm
gc = 32.2;%(lbm*ft/s^2)/lbf

%% Freestream
%M0 and Alt are given as flight conditions (unitless, feet)
[T0, a0, P0, ~] = atmosisa_EE(Alt); %degR, ft/s, psia
[~, trinv, PIrinv, ~, ~] = flowisentropic(gamma, M0); % Ram ratios, unitless

tr = 1/trinv;
PIr = 1/PIrinv;

Tt0 =tr*T0; % degR
Pt0 = PIr*P0; %psia
V0 = a0*M0; %ft/s
%% Diffuser
td = 1;
PId = 1;

Tt2 = td*Tt0; % Determined via enegy balance (ht0 = ht2)
Pt2 = PId*Pt0; %td^(gamma/(gamma - 1))

%mdotc = mdot0/(1+alpha);

%mdotf = mdot0 - mdotc;
%% Fan
%PIf is given as a design choice
Pt13 = Pt2*PIf;
tf = PIf^((gamma-1)/gamma);
Tt13 = Tt2*tf;

%% Compressor
% PIc is given as given as a design choice
tc = PIc^((gamma-1)/gamma); %Total temp ratio over compressor

Tt3 = tc*Tt2; %degR, total temp after compressor
Pt3 = PIc * Pt2; %psia, total pres after compressor

%% Burner 
%Tt4 is given

tb = Tt4/Tt3; % Total temp ratio across burner
tl = Tt4/T0; % Ratio between total Temp after burner and static freestream
PIb = 1; %Isobaric heat addition

Pt4 = Pt3*PIb; %Total pres after burner, degR
f = cp/Hpr*(Tt4 - Tt3); %mdotfuel/mdotc
%mdotfuel = f*mdotc;
%% Turbine
tt = 1-tr/tl*(tc-1+alpha*(tf-1));%Power blance from compressor and fan
Tt5 = Tt4*tt;
PIt = tt^(gamma/(gamma-1)); 
Pt5 = PIt * Tt4;

%% Core Nozzle
tn = 1; %Adiabatic, No work;
PIn = 1;

Tt9 = tn * Tt5; %degR, Exit total Temp
Pt9 = PIn * Pt5; %psia, Exit total Pressure
P9 = P0; %Subsonic/Ideal Expansion
%% Fan Nozzle
tfn = 1; %No losses
PIfn = 1;
Tt19 = Tt13*tfn;
Pt19 = Pt13*PIfn;
P19 = P0; %Ideal/subsonic nozzle
%% Fan Exit Conditions
P19_Pt19 = P19/Pt19; % Also = 1/(PIf*PIr)
if ((P19_Pt19 < 1) && (P19_Pt19 > 0))
    [M19, T19_Tt19, ~, ~, ~] = flowisentropic(gamma, P19_Pt19, 'pres');
else
    M19 = NaN; 
    T19_Tt19 = NaN;
end
T19 = T19_Tt19 * Tt19;
a19 = sqrt(gamma*R*T19*gc);
V19 = a19*M19; % exit velocity

%% Core Exit Conditions
P9_Pt9 = P9/Pt9; %Static /Total temp ratio
if ((P9_Pt9 < 1) && (P9_Pt9 > 0) && isreal(P9_Pt9))
    M9 = sqrt(2/(gamma-1)*(tt*tc*tr-1));
    T9 = T0*tl/(tr*tc);
else
    M9 = NaN; 
    T9 = NaN;
end
    

a9 = sqrt(gamma*R*T9*gc);
V9 = a9*M9; % exit velocity
%% Efficiencies

nt = 1-1/(tr*tc);%(alpha*V19^2+V9^2-alpha*(alpha+1)*V0^2)/(f*Hpr); %thermal efficiency
np = 2*V0*(alpha*(V19-V0)+(V9-V0))/(alpha*(V19^2-V0^2)+(V9^2-V0^2)); %propulsive efficiency
no = np * nt; %overall efficiency
F_mdot0 = 1/gc*1/(alpha+1)*((V9-V0)+alpha*(V19-V0));
S = f/((alpha+1)*F_mdot0)*3600;

%% Output
Pmat(1) = P0;
Pmat(6) = P9;
Pmat(8) = P19;

Tmat(1) = T0;
Tmat(6) = T9;
Tmat(8) = T19;

Mmat(1) = M0;
Mmat(6) = M9;
Mmat(8) = M19;

Vmat(1) = V0;
Vmat(6) = V9;
Vmat(8) = V19;

Emat = [nt np no F_mdot0 S];

Ptmat = [Pt0 Pt2 Pt3 Pt4 Pt5 Pt9 Pt13 Pt19];
Ttmat = [Tt0 Tt2 Tt3 Tt4 Tt5 Tt9 Tt13 Tt13];
end

