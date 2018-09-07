function [F, S] = uninstFS(mdot_0, mdot_fan, mdot_9, mdot_fuel, M_0, Alt, V_fan, V_9, unit)
    if unit == 'EE'
        gc = 32.2;
        [T, a, P, rho] = atmosisa(Alt);
    elseif unit == 'SI'
        gc = 1;
        Alt = convlength(Alt, 'ft', 'm');
        [T, a, P, rho] = atmosisa(Alt);
        a = convlength(a, 'm', 'ft');
    else
        return
    end
    V_0 = M_0*a;
    F = ((mdot_9*V_9 + mdot_fan*V_fan)-mdot_0*V_0)/gc;
    S = mdot_fuel/F;
end