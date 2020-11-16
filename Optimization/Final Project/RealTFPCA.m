function [Tmat,Pmat, Ttmat, Ptmat, Mmat, Vmat, Emat, PImat, tmat] = RealTFPCA(M0,Alt, Tt4, PIf, alpha, PIc,PIdmax,ef, ec,nb, nm,PIb,et,PIn, PIfn)
%[Tmat,Pmat, Ttmat, Ptmat, Mmat, Vmat, Emat] = RealTFPCA(M0,Alt, Tt4, PIf, alpha, PIc,PIdmax,ec,nb, nm,PIb,et,PIn, PIfn)
%[degR, psia, degR,psia, unitless, ft/s]=fcn(unitless, ft, degR,3 unitless design inputs, 6 unitless figs of merit))
%% Constants
Rc = 53.35; %(lbf ft)/(lbm degR)
Rt = 52.96; %(lbf ft)/(lbm degR)
cpc = 0.24; % Btu/(lbm ft/s^2)
cpt = .295; % Btu/(lbm ft/s^2)
gammac = 1.4; %unitless
gammat = 1.3; %unitless
Hpr = 18400; %Btu/lbm
gc = 32.2;%(lbm*ft/s^2)/lbf

%% Freestream
%M0 and Al  t are given as flight conditions (unitless, feet)
[T0, a0, P0, ~] = atmosisa_EE(Alt); %degR, ft/s, psia

tr = 1+(gammac-1)/2*M0^2;
PIr = tr^(gammac/(gammac-1));

Tt0 =tr*T0; % degR
Pt0 = PIr*P0; %psia
V0 = a0*M0; %ft/s

%% Diffuser

if M0 > 1
    nr = 1-0.075*(M0-1)^(1.35);
else
    nr = 1;
end

PId = PIdmax*nr;
td = PId^((gammac-1)/(gammac));

Tt2 = td*Tt0; %Figures of merit
Pt2 = PId*Pt0; %Figures of merit
%% Fan
Pt13 = PIf*Pt2;
tf = PIf^((gammac-1)/(gammac*ef));
Tt13 = Tt2*tf;
nf = (PIf^((gammac-1)/(gammac))-1)/(tf-1);

%% Compressor
% PIc is given as given as a design choice
tc = PIc^((gammac-1)/(gammac*ec)); %Total temp ratio over compressor (with fig of merit)
nc = (PIc^((gammac-1)/(gammac))-1)/(tc-1);
Tt3 = tc*Tt2; %degR, total temp after compressor
Pt3 = PIc * Pt2; %psia, total pres after compressor

%% Burner 
%Tt4 and PIb are given

tb = Tt4/Tt3; % Total temp ratio across burner
tl = cpt*Tt4/(cpc*T0); % Ratio between total Temp after burner and static freestream
ht3 = cpc*Tt3;

f = fIterate(ht3,nb,Hpr,Tt4);

Pt4 = Pt3*PIb; %Total pres after burner, degR

%% Turbine
tt = 1-1/(nm*(1+f))*tr/tl*(tc-1+alpha*(tf-1));
Tt5 = Tt4*tt;%degR
PIt = tt^(gammat/((gammat-1)*et));% 
nt = (1-tt/(1-tt^(1/et)));
Pt5 = PIt * Tt4;

%% Core Nozzle
P9 = P0/((((gammat+1)/2)^(gammat/(gammat-1)))/(PIr*PId*PIc*PIb*PIt*PIn));
Pt9 = P9*((gammat+1)/2)^(gammat/(gammat-1));
if P9 < P0 %If not choked
    P9 = P0;
    Pt9 = P0*PIr*PId*PIc*PIb*PIt*PIn;
end

M9 = sqrt(2/(gammat-1)*((Pt9/P9)^((gammat-1)/gammat)-1));

Tt9 = T0*tl*tt*(cpc/cpt);

T9 = T0*tl*tt*(cpc/cpt)/((Pt9/P9)^((gammat-1)/gammat));

V9 = a0*M9*sqrt((gammat*Rt*T9)/(gammac*Rc*T0));

%% Fan Nozzle
P19 = P0/((((gammac+1)/2)^(gammac/(gammac-1)))/(PIr*PId*PIf*PIfn));
Pt19 = P19*((gammac+1)/2)^(gammac/(gammac-1));
if P19 < P0 %If not choked
    P19 = P0;
    Pt19 = P0*PIr*PId*PIf*PIfn;
end

M19 = sqrt(2/(gammac-1)*((Pt19/P19)^((gammac-1)/gammac)-1));

Tt19 = T0*tr*tf;

T19 = T0*tr*tf/((Pt19/P19)^((gammac-1)/gammac));

V19 = a0*M19*sqrt((T19)/(T0));

%% Performance
F_mdot0 = 1/(1+alpha)*(a0/gc*((1+f)*V9/a0-M0+(1+f)*Rt/Rc*(T9/T0)/(V9/a0)*(1-P0/P9)/gammac))+alpha/(1+alpha)*a0/gc*(V19/a0-M0+(T19/T0)/(V19/a0)*(1-P0/P19)/gammac);

nt = a0^2*((1+f)*(V9/a0)^2+alpha*(V19/a0)^2-(1+alpha)*M0^2)/(2*gc*f*Hpr*778); %thermal efficiency

np = 2*gc*V0*(1+alpha)*F_mdot0/(a0^2*((1+f)*(V9/a0)^2+alpha*(V19/a0)^2-(1+alpha)*M0^2)); %propulsive efficiency

no = np * nt; %overall efficiency

S = f/((1+alpha)*F_mdot0)*3600;

FR = ((1+f)*V9/a0-M0+(1+f)*Rt/Rc*(T9/T0)/(V9/a0)*(1-P0/P9)/gammac)/(V19/a0-M0+(T19/T0)/(V19/a0)*(1-P0/P19)/gammac);
%% Function Output
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

Emat = [nt np no F_mdot0 S nc nf nt];

PImat = [PId PIt PIr];

tmat = [tf tc tt tl tr];

Ptmat = [Pt0 Pt2 Pt3 Pt4 Pt5 Pt9 Pt13 Pt19];
Ttmat = [Tt0 Tt2 Tt3 Tt4 Tt5 Tt9 Tt13 Tt19];

