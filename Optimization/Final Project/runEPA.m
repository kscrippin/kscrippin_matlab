clear 
close all
clc

M0 = 0.1;
Alt = 5000;
Tt4 = 3200;

mdot0r = 100; %lbm/s
M0r = 0.8;
Altr = 36000; % ft
Tt4r = 3200; %degR
PIfr = 3;
PIdmax = 1;
alphar = 7.6;
PIcr = 12;
ef = 1;
ec = 1;
nb = 1;
PIb = 1;
et = 1;
PIn = 1;
PIfn = 1;
nm = 1;

[Tmatr,Pmatr, Ttmatr, Ptmatr, Mmatr, Vmatr, Ematr, PImatr, tmatr] = RealTFPCA(M0r,Altr, Tt4r, PIfr, alphar, PIcr,PIdmax,ef, ec,nb, nm,PIb,et,PIn, PIfn);

nc = Ematr(6);
nf = Ematr(7);
nt = Ematr(8);
PIdr = PImatr(1);
PItr = PImatr(2);
PIrr = PImatr(3);
tfr = tmatr(1);
tcr = tmatr(2);
ttr  = tmatr(3);
tlr = tmatr(4);
trr = tmatr(5);
M9r = Mmatr(6);
M19r = Mmatr(8);

[F,S,no] = TFEPA(M0,Alt,Tt4,PIdr,PItr,PIrr,tfr,tcr,ttr,tlr,trr,M0r,M9r,M19r,Altr, Tt4r, PIfr, alphar, PIcr,PIdmax,ef, ec,nb, nm,PIb,et,PIn,PIfn,nc,nf,nt, mdot0r)