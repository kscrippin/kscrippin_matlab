function [Temp, Va, Pres, density] = atmosisa_EE(height)
%[Temp, Va, Pres, density] = atmosisa_EE(height)
%height is in ft, Temp degR, Pres psi, a ft/s, density lbm/ft^3

height_SI = convlength(height, 'ft', 'm');

[T, a, P, rho] = atmosisa(height_SI);

Temp = convtemp(T, 'K', 'R');
Pres = convpres(P, 'Pa', 'psi');
Va = convvel(a, 'm/s', 'ft/s');
density = convdensity(rho, 'kg/m^3', 'lbm/ft^3');