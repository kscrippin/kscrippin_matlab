clear;
clc;

%-------------User Input------------
mdot_0 = [472];     %lbm/s or kg/s
mdot_fan = [412];   %lbm/s or kg/s
mdot_fuel = [1.028]; %lbm/s or kg/s
mdot_9 = [(mdot_0-mdot_fan)]; %lbm/s or kg/s
M_0 = [.82]; %unitless
Alt = [32000]; %ft or m
V_9 = [1004]; %ft/s or m/s
V_fan = [1298]; %ft/s or m/s
unit = ['EE']; %ft/s or m/s

%------------------------------------
m = length(mdot_0); %please ensure all input arrays above are of the same length
F = zeros(1, m);
S = zeros(1, m);
for n = (1:m)
    [F(n), S(n)] = uninstFS(mdot_0(n), mdot_fan(n), mdot_9(n), mdot_fuel(n), M_0(n), Alt(n), V_fan(n), V_9(n), unit(n));
end