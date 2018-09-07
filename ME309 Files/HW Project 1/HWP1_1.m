%Script to run multiple iterations of the aspeed funtion

%-------------User Input------------------
temperature = [405 3240 225 1800]; % degR OR K
gasConstant = [53.36 52.95 287 285]; %(Ft*lbf/lbm*degR) OR (J/kg*K)
heatCapacityRatio = [1.4 1.3 1.4 1.3]; %Unitless
units = ['EE' 'EE' 'SI' 'SI']; %'EE' OR 'SI'
%----------------------------------------

len = length(temperature);
acousticSpeed = ones(1, len);
Est_Alt = ones(1, len);

for n = (1:len)
   [acousticSpeed(n), Est_Alt(n)] = aspeed(temperature(n), gasConstant(n), heatCapacityRatio(n), units(2*n));
end

disp(acousticSpeed);