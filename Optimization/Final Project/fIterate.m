function f = fIterate(ht3,nb,Hpr,Tt4)
%Finds fuel air ratio for Real PCA
%   [f,ht4] = fIterate(ht3,nb,Hpr,ht4guess,fguess)
%--------------------
convcrit = 0.0001;
%--------------------
s = open('fuelair.mat');

cpt = 1.3;
ht4guess = cpt*Tt4;
fguess = 0.0169; %mdotf/mdot0
flag = true;
Temptable = s.Temptable;
ht4table = s.ht4table;
ftable = s.ftable;

while flag
    
    n =2;
    while (Temptable(n,1) <= Tt4)
        n = n+1;
    end
    m = 2;
    while ftable(1,m) <= fguess
        m = m+1;
    end
    
    ht4check = ((fguess-ftable(m-1))/(ftable(m)-ftable(m-1)))*(ht4table(n,m)-ht4table(n,m-1))+ht4table(n,m-1);
    fcheck = (ht4check-ht3)/(nb*Hpr-ht4check);
    if (fcheck - fguess) < convcrit
        flag = false;
    end
        fguess = fcheck;
end
f = fcheck;
ht4 = ht4check;
end

