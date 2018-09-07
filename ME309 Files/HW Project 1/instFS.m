function [T, TSFC] = instFS(mdot_0, mdot_fan, mdot_9, mdot_fuel, M_0, Alt, V_fan, V_9, unit)
    %Reports Thrust and TSFC graphically and in an array
    %[T, TSFC] = instFS(mdot_0, mdot_fan, mdot_9, mdot_fuel, M_0, Alt, V_fan, V_9, unit)
    %mdot - lbm/s or kg/s, Alt - ft or m, V - ft/s or m/s, unit - 'EE' 'SI'
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
    S = mdot_fuel/F * 3600;
    
    phi_in = linspace(0,0.1);
    phi_noz = linspace(0, .05);
    
    T = F.*(1 - phi_in - phi_noz);
    TSFC = S./(1 - phi_in - phi_noz);
    perPhi = linspace(0, 100); %Percentage of drag coefficients
    
    figure
    plot(perPhi, T, '-k'); 
    title("Thrust vs Drag Coefficient (\phi_in = 0.1, \phi_nox = 0.05)")
    xlabel("Percentage of Total Drag Ceofficient")
    ylabel("Installed Thrust (lb_f)")
    
    figure
    plot(perPhi, TSFC, '-b');
    title("TSFC vs Drag Coefficient (\phi_in = 0.1, \phi_nox = 0.05)")
    xlabel("Percentage of Total Drag Ceofficient")
    ylabel("TSFC (lb_m/(hr*lb_f))")
    
    figure
    subplot(2,1,1);
    plot(perPhi, T, '-k');
    title("Thrust and TSFC vs Drag Coefficient (\phi_i = 0.1, \phi_n = 0.05)")
    xlabel("Percentage of Total Drag Ceofficient")
    ylabel("Installed Thrust (lb_f)")
    subplot(2,1,2);
    plot(perPhi, TSFC, '-b');
    xlabel("Percentage of Total Drag Ceofficient")
    ylabel("TSFC (lb_m/(hr*lb_f))")
end