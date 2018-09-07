function [T, TSFC] = instFS(mdot_0, mdot_fan, mdot_9, mdot_fuel, M_0, Alt, V_fan, V_9, unit)
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
    phi_in = linspace(0,0.1);
    phi_noz = linspace(0, .05);
    T = F.*(1 - phi_in - phi_noz);
    sArray = S.*ones(1, length(phi_in));
    TSFC = S./(1 - phi_in - phi_noz);
    perPhi = linspace(0, 100);
    figure
    plot(perPhi, T, '-k');
    figure
    plot(perPhi, TSFC, '-b');
    figure
    subplot(2,1,1);
    plot(perPhi, T, '-k');
    subplot(2,1,2);
    plot(perPhi, TSFC, '-b');
end