clear 
close all
clc

M0 = 0.8;
Alt = 36000; % ft
Tt4 = 3200; %degR
PIf = 3;
PIdmax = 1;
alpha = 7.6;
PIc = 12;
ef = 1;
ec = 1;
nb = 1;
PIb = 1;
et = 1;
PIn = 1;
PIfn = 1;
nm = 1;

[Tmat,Pmat, Ttmat, Ptmat, Mmat, Vmat, Emat] = RealTFPCA(M0,Alt, Tt4, PIf, alpha, PIc,PIdmax,ef, ec,nb, nm,PIb,et,PIn, PIfn);
