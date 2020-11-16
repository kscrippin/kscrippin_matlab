function [F,S,no] = TFEPA(M0,Alt,Tt4,PIdr,PItr,PIrr,tfr,tcr,ttr,tlr,trr,M0r,M9r,M19r,Altr, Tt4r, PIfr, alphar, PIcr,PIdmax,ef, ec,nb, nm,PIb,et,PIn,PIfn,nc,nf,nt,mdot0r)
%% TurboFan Engine Performance Analysis
%This code takes flight condition inputs based on Ideal PCA (soon to be Real) and generates a EPA based on the altered Mach number, altitude, and throttle setting. (For an ideal, single spool, turbofan engine.)

Rc = 53.35; %(lbf ft)/(lbm degR)
Rt = 52.96; %(lbf ft)/(lbm degR)
cpc = 0.24; % Btu/(lbm ft/s^2)
cpt = .295; % Btu/(lbm ft/s^2)
gammac = 1.4; %unitless
gammat = 1.3; %unitless
Hpr = 18400; %Btu/lbm
gc = 32.2;%(lbm*ft/s^2)/lbf
    
    % Choices
        %Flight Paramters

        [T0, a0, P0, ~] = atmosisa_EE(Alt); %[deg R, ft/s, psi, lbm/ft^3]
        [T0r, a0r, P0r, ~] = atmosisa_EE(Altr);

        %Gas Properties
        gamma_c = 1.4;
        gamma_t = 1.3;

        c_pc = 0.24; %Btu/(lbm*degR)
        c_pt = 0.295; %Btu/(lbm*degR)

        Rc = 53.35;
        Rt = 52.96;

        %Fuel
        h_PR = 18400; 


    %% conditions


    tr = 1+(gammac-1)/2*M0^2;
    PIr = tr^(gammac/(gammac-1));

    V0 = a0*M0; %ft/s

    if M0 > 1
        nr = 1-0.075*(M0-1)^(1.35);
    else
        nr = 1;
    end

    PId = PIdmax*nr;
    td = PId^((gammac-1)/(gammac));
    tl = cpt*Tt4/(cpc*T0);

    ttcheck = ttr;
    PItcheck = PItr;
    tfcheck = tfr;

    tt = 0;

    while abs(tt-ttcheck) > 0.0001

    tt = ttcheck;
    PIt = PItcheck;
    tf = tfcheck;


    tc = 1+(tl/tr)/(tlr/trr)*tfr/tf*(tcr-1);

    PIc = (1+(tc-1)*nc)^(gammac/(gammac-1));

    PIf = (1+(tf-1)*nf)^(gammac/(gammac-1));

    Pt19 = P0*PIr*PId*PIf*PIfn;

    if (Pt19/P0) < (((gammac+1)/2)^(gammac/(gammac-1)))
        P19 = P0;
    else
        P19 = Pt19/(((gammac+1)/2)^(gammac/(gammac-1)));
    end
    M19 = sqrt(2/(gammac-1)*((Pt19/P19)^((gammac-1)/gammac)-1));

    Pt9 = P0*PIr*PId*PIf*PIc*PIb*PIt*PIn;

    if (Pt9/P0) < (((gammat+1)/2)^(gammat/(gammat-1)))
        P9 = P0;
    else
        P9 = Pt9/(((gammat+1)/2)^(gammat/(gammat-1)));
    end

    M9 = sqrt(2/(gammat-1)*((Pt9/P9)^((gammat-1)/gammat)-1));

    MFP9 = sqrt(gammat*gc/Rt)*M9/(1+(gammat-1)/2*M9^2)^((gammat+1)/(2*(gammat-1)));
    MFP9r = sqrt(gammat*gc/Rt)*M9r/(1+(gammat-1)/2*M9r^2)^((gammat+1)/(2*(gammat-1)));
    
    MFP19 = sqrt(gammac*gc/Rc)*M19/(1+(gammac-1)/2*M19^2)^((gammac+1)/(2*(gammac-1)));
    MFP19r = sqrt(gammac*gc/Rc)*M19r/(1+(gammac-1)/2*M19r^2)^((gammac+1)/(2*(gammac-1)));

    alpha = alphar*(PIcr/PIc*sqrt((tl*(tr*tf))/(tlr/(tr*tf))))*MFP19r/MFP19;

    tfcheck = 1+(1-tt)/(1-ttr)*(tl/tr)/(tlr/trr)*(1+alphar)/(1+alpha)*(tfr-1);

    ttcheck = 1-nt*(1-PIt^((gammat-1)/gammat));

    PItcheck = PItr*sqrt(tt/ttr)*MFP9r/MFP9;
    end

    mdot0 = mdot0r*(1+alpha)/(1+alphar)*(P0*PIr*PId*PIf*PIc)/(P0r*PIrr*PIdr*PIfr*PIcr)*sqrt(Tt4r/Tt4)

    Tt3 = T0*tr*td*tf*tc;

    ht3 = cpc*Tt3;

    f = fIterate(ht3,nb,Hpr,Tt4);

    T9 = T0 * tl*tt*(cpc/cpt)/((Pt9/P9)^((gammat-1)/gammat));
    T19 = T0 * tr*tf*(cpc/cpt)/((Pt19/P19)^((gammac-1)/gammac));

    V9 = a0*M9*sqrt((gammat*Rt*T9)/(gammac*Rc*T0));
    V19 = a0*M19*sqrt((gammat*Rt*T19)/(gammac*Rc*T0));

    F_mdot0 = 1/(1+alpha)*(a0/gc*((1+f)*V9/a0-M0+(1+f)*Rt/Rc*(T9/T0)/(V9/a0)*(1-P0/P9)/gammac))+alpha/(1+alpha)*a0/gc*(V19/a0-M0+(T19/T0)/(V19/a0)*(1-P0/P19)/gammac);
    F = F_mdot0*mdot0;

    S = f/((1+alpha)*(F_mdot0));
    TSFC = S*3600;

    nt = a0^2*((1+f)*(V9/a0)^2+alpha*(V19/a0)^2-(1+alpha)*M0^2)/(2*gc*f*Hpr*778) %thermal efficiency

    np = 2*gc*V0*(1+alpha)*F_mdot0/(a0^2*((1+f)*(V9/a0)^2+alpha*(V19/a0)^2-(1+alpha)*M0^2)) %propulsive efficiency

    no = np * nt; %overall efficiency
  end