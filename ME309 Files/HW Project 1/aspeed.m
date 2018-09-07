function a = aspeed(T, R, Gamma, Unit)
    %acoustic_speed = aspeed(Temperature (absolute), Gas Constant, Heat
    %Capacity Ratio, Unit System ('EE' or 'SI')
    if Unit == 'EE'
        Gc = 32.2;
    elseif Unit == 'SI'
        Gc = 1;
    else
        disp('Error, please specifiy if the unit system is SI or EE in a string array');
        a = NaN;
        return
    end
    a = sqrt(Gc*Gamma*R*T);
end