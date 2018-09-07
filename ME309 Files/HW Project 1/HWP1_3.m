%-------------User Input------------
mdot_0 = 472;     %lbm/s or kg/s
mdot_fan = 412;   %lbm/s or kg/s
mdot_fuel = 1.028; %lbm/s or kg/s
mdot_9 = (mdot_0-mdot_fan); %lbm/s or kg/s
M_0 = .82; %unitless
Alt = 32000; %ft or m
V_9 = 1004; %ft/s or m/s
V_fan = 1298; %ft/s or m/s
unit = 'EE'; %ft/s or m/s
%------------------------------------

[T, TSFC] = instFS(mdot_0, mdot_fan, mdot_9, mdot_fuel, M_0, Alt, V_fan, V_9, unit);
% lbf, lbm/(hr*lbf)