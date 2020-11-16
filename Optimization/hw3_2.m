%% Homework 3 Problem 2

clear
close all
clc

% syms s1 s2 s3 x1 x2 L1 L2
% 
% s = [s1 s2;s2 s3];
% a = [0 1;0 0];
% b = [0;1];
% r = 2;
% q = [1 0;0 1];
% x = [x1; x2];
% L = [L1; L2];
% 
% sdot = -1*(s*a+a'*s-s*b*r^-1*b'*s+q)
% u = -1*r^-1*b'*L;
% 
% k = -1*r^-1*b'*s
% 
% sdotb = s*a+a'*s-s*b*r^-1*b'*s+q == 0;
% 
% S = solve([sdotb(1,1),sdotb(2,1),sdot(2,2)],[s1,s2,s3]);

[T, S] = ode45(@(t,s)eom(t,s),flip(linspace(0,10,1001)),[1,0,1]);
S = flip(S);
T = flip(T);

[T2, S2] = ode45(@(t,s)eom2(t,s,T,S),linspace(0,10,1001),[2,0]);
[T3, S3] = ode45(@(t,s)eom3(t,s),linspace(0,10,1001),[2,0]);

figure(1)
plot(T,S(:,1),'-k',T,S(:,2),'-b',T,S(:,3),'-r')
title('S(t)')
xlabel('time, s')
ylabel('S')
legend('S1','S2','S3')

figure(2)
plot(T2,S2(:,1),'-k',T2,S2(:,2),'-b',T3,S3(:,1),'--r',T3,S3(:,2),'--g')
title('X(t)')
xlabel('time, s')
ylabel('X')
legend('X1','X2','X1 Alg','X2 Alg')

figure(3)
U = -0.5.*(S(:,2).*S2(:,1)+S(:,3).*S2(:,2));

Ualg = -0.5.*(1.414.*S2(:,1)+2.767.*S2(:,2));

plot(T2,U,'-b',T3,Ualg,'--r');
title('U(t)')
xlabel('time, s')
ylabel('U')
legend('U1','U1 Alg')

function ds = eom(~,s)
% s(1) = S1, s(2) = S2, s(3) = S3
ds(1) = 0.5*s(2)^2-1;
ds(2) = s(2)*s(3)/2-s(1);
ds(3) = 0.5*s(3)^2-2*s(2)-1;
ds = ds';
end

function ds = eom2(t,s,T,S)
% s(1) = x1, s(2) = x2
s2 = interp1(T,S(:,2),t);
s3 = interp1(T,S(:,3),t);

ds(1) = s(2);
ds(2) = -s2/2*s(1)+-s3/2*s(2);

ds = ds';
end

function ds = eom3(t,s)
% s(1) = x1, s(2) = x2
s2 = 1.414;
s3 = 2.767;

ds(1) = s(2);
ds(2) = -s2/2*s(1)+-s3/2*s(2);

ds = ds';
end