function [a, Alt_est] = aspeed(T, R, Gamma, Unit)
    %acoustic_speed = aspeed(Temperature (absolute), Gas Constant, Heat
    %Capacity Ratio, Unit System ('EE' or 'SI')
    if Unit == 'E'
        Gc = 32.2;
        a_isa = 1116;
    elseif Unit == 'I'
        Gc = 1;
        a_isa = 340.2941;
    else
        disp('Error, please specifiy if the unit system is SI or EE in a string array');
        a = NaN;
        return
    end
    
    a = sqrt(Gc*Gamma*R*T);
    Alt_est = 0;
   
   if Gamma == 1.4
        while a_isa > a;
            if Unit == 'E'
                [Temp, a_isa, P, rho] = atmosisa_EE(Alt_est);
                Alt_est = Alt_est + 1000;
            elseif Unit == 'I' && Gamma == 1.4
                [Temp, a_isa, P, rho] = atmosisa(Alt_est);
                Alt_est = Alt_est + 1000;
            else
            Alt_est = NaN;
            end
        end
   else
        Alt_est = NaN;
   end
end