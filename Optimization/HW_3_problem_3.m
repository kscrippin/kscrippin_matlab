clear
close all
clc
%% Part B
% solve for k and eigen vector values 
A =-1*[0 -1 1 0
       -1 1 0 0
        1 0 0 1
       0 1 1 -1];
% v is "k" and s is eigen vector

[V, s] = eig(A)
% solving for C coefficients

syms A B C D t
% this is the "c" coeffiecient expression
x = A.*V(:,1)*exp(s(1,1).*t)+B.*V(:,2).*exp(s(2,2).*t)+C.*V(:,3).*exp(s(3,3).*t)+D.*V(:,4).*exp(s(4,4).*t);

%we then sub t as 0 in all the x functions below to indicate initial time
x_0 = subs(x,t,0);
eqn(1) = 0 == x_0(1);
eqn(2) = 0 == x_0(2);

% x 2 has two forms, one at t = 0 and one at t = tf which is equation 3
%equation 4 is lambd two which is defined by x1 at x1 at final time
x_f = simplify(subs(x,t,2));
eqn(3) = 1 == x_f(2);
eqn(4) = 10*x_f(1) == x_f(3);
% solving for the coeficients...
k = [A B C D];
ks = solve(eqn,k);
kd(1) = double(ks.A);
kd(2) = double(ks.B);
kd(3) = double(ks.C);
kd(4) = double(ks.D);
kd

%sub the coefficients back into the original x equation to create a scalar
x = subs(x,k,kd);

%integrate the 1/2 x^2 + x^2 +u^2 where u = -lambda1
int_H = int((x(1))^2+(x(2))^2+(-1*x(3))^2,t,0,2);

% solve for nu and cost
tf = 2;
nu = (kd(1).*(V(4,1)*exp(s(1,1).*tf))) + (kd(2).*(V(4,2)*exp(s(2,2).*tf))) ...
    + (kd(3).*(V(4,3)*exp(s(3,3).*tf))) + (kd(4).*(V(4,4)*exp(s(4,4).*tf)));

x1_tf = (kd(1).*(V(1,1)*exp(s(1,1).*tf))) + (kd(2).*(V(1,2)*exp(s(2,2).*tf))) ...
    + (kd(3).*(V(1,3)*exp(s(3,3).*tf))) + (kd(4).*(V(1,4)*exp(s(4,4).*tf)));

J = double(5*((x1_tf)^2) + 0.5*(int_H));

% important t information, gotta bring down into the bv4c function for
% initl guess
t_1 = linspace(0,2,50);
% solve for each individual x1, x2, x3, x4 as an array of values...

x1 = (kd(1).*(V(1,1)*exp(s(1,1).*t_1))) + (kd(2).*(V(1,2)*exp(s(2,2).*t_1))) ...
    + (kd(3).*(V(1,3)*exp(s(3,3).*t_1))) + (kd(4).*(V(1,4)*exp(s(4,4).*t_1)));

x2 = (kd(1).*(V(2,1)*exp(s(1,1).*t_1))) + (kd(2).*(V(2,2)*exp(s(2,2).*t_1))) ...
    + (kd(3).*(V(2,3)*exp(s(3,3).*t_1))) + (kd(4).*(V(2,4)*exp(s(4,4).*t_1)));

x3 = (kd(1).*(V(3,1)*exp(s(1,1).*t_1))) + (kd(2).*(V(3,2)*exp(s(2,2).*t_1))) ...
    + (kd(3).*(V(3,3)*exp(s(3,3).*t_1))) + (kd(4).*(V(3,4)*exp(s(4,4).*t_1)));

x4 = (kd(1).*(V(4,1)*exp(s(1,1).*t_1))) + (kd(2).*(V(4,2)*exp(s(2,2).*t_1))) ...
    + (kd(3).*(V(4,3)*exp(s(3,3).*t_1))) + (kd(4).*(V(4,4)*exp(s(4,4).*t_1)));

u = -1 * ((kd(1).*(V(3,1)*exp(s(1,1).*t_1))) + (kd(2).*(V(3,2)*exp(s(2,2).*t_1))) ...
    + (kd(3).*(V(3,3)*exp(s(3,3).*t_1))) + (kd(4).*(V(3,4)*exp(s(4,4).*t_1))));

%% ODE45 part
% solve for the initial x1 and x2 with ode 45
[Tode,Xode] = ode45(@(tode,xode)eom(tode,xode),[0 2],[.1,.12]);
u_guess = linspace(.1,.1,41);
%solve for the initial lambda1 and lambda2 with ode 45
[T2,Lambda] = ode45(@(tode,xode)LagrangeEom(tode,xode,Tode,Xode),Tode(end:-1:1),[-1,3]);


%% BVP4C Section Part d and Part h
% initial guesses need to be checked , solinit. x should be time, but i am
% unsure whether we can use the calculated x values for solinit.y


% Possibly WRONG *****************************************************************************************************
solx = (Tode');
soly = ([Xode(:,1), Xode(:,2),Lambda(:,1),Lambda(:,2)])';
solinit.x = solx;
solinit.y = soly;

%BVP4C calling statement:

sol = bvp4c(@(te,se)ELeom(te,se),@bcfunc,solinit);

%BVP4C Free time calling statement:
sol2 = bvp4c(@(te,se)ELeom(te,se),@bcfunc2,solinit);

%solve for optimal cost using trapz for integral part for part h
u2 = -1.*(sol2.y(3,:));
%call ode 45 to integrate (x1^2 + x2^2 + u^2)
t41 = linspace(0,2,41)';
% [J2] = ode45(@(utode,uode)eom2(sol2.x,sol2.y(1,:),sol2.y(2,:),sol2.y(3,:),Tode,Xode),u2(end:-1:1),t41');
Y = (sol2.y(1,:).^2) + (sol2.y(2,:).^2) + (sol2.y(3,:).^2);
J2 = 5.*((sol2.y(1,end).^2)) + 0.5.*(trapz(t41,Y));

%% Sub Plots
% Part C
figure(1)
subplot(2,1,1)
plot(t_1,x1,'-k',t_1,x2, '-b')
xlabel('time, sec')
ylabel('x values')
legend('x_1', 'x_2','Location','Northwest')
title(" Optimal cost: "+num2str(J))
subplot(2,1,2)
plot(t_1,u,'--r')
xlabel('time, sec')
ylabel('u values')

%part e

figure(2)
subplot(2,1,1)
plot(solinit.x,solinit.y(1,:),'-k',solinit.x,solinit.y(2,:),'-b',sol.x,sol.y(1,:),'-r',sol.x,sol.y(2,:),'-g')
xlabel('time, sec')
ylabel('x values')
legend('x_1 initial guess', 'x_2 initial guess','x_1 bvp4c','x_2 bvp4c','Location','Northwest')
subplot(2,1,2)
plot(solinit.x,u_guess,'-r',sol.x,(-1.*sol.y(3,:)),'-g')
xlabel('time, sec')
ylabel('x values')
legend('u initial guess', 'u bvp4c','Location','Northwest')

%part f
figure(3)
subplot(2,1,1)
% analytical solution
plot(t_1,x1,'ok',t_1,x2, 'ob',sol.x,sol.y(1,:),'-r',sol.x,sol.y(2,:),'-g')
xlabel('time, sec')
ylabel('x values')
legend('x_1 analytical', 'x_2 analytical','x_1 bvp4c', 'x_2 bvp4c','Location','Northwest')
% bvp4c solution
subplot(2,1,2)
plot(t_1,(-1.*x3),'ok',sol.x,(-1.*sol.y(3,:)),'-g')
xlabel('time, sec')
ylabel('x values')
legend('u analytical', 'u bvp4c','Location','Northwest')

% part h
figure(4)
subplot(2,1,1)
plot(sol2.x,sol2.y(1,:),'-k',sol2.x,sol2.y(2,:),'-b')
xlabel('time, sec')
ylabel('x values')
legend('x_1', 'x_2','Location','Northwest')
title(" Optimal cost: "+num2str(J2))
%
subplot(2,1,2)
plot(sol2.x,u2,'--r')
xlabel('time, sec')
ylabel('u values')

%% Functions
function ds = eom(tode,sode)
ds(1) = -1*(sode(2) - 0.1);
ds(2) = -1*(sode(1) - sode(2));
ds =ds';
end

function J2 = eom2(utode, x1,x2,x3,T,UX)
XY = interp1(T,UX,utode);
J2 = double(5*((x1).^2) + 0.5*(x1.^2 + x2.^2 + x3.^2));
J2 =J2';
end

function ds = LagrangeEom(t,s,T,SXY)
% s(1) = lambda_1, s(2) = lambda_2
XY = interp1(T,SXY,t);
ds(1) = -1*(XY(1) + s(2));
ds(2) = -1*(XY(2) + s(1) - s(2));
ds = ds';
end

function ds = ELeom(~,se)
% s(1) = x_1, s(2) = x_2, s(3) = Lambda_1, s(4) = Lambda_2
%u = -1*se(3); %Optimal Control Law
ds(1) = se(2)-se(3); %xdot1
ds(2) = se(1)-se(2); %xdot2
ds(3) = -1*(se(1)+se(4)); %lambdadot1
ds(4) = -1*(se(2)+se(3)-se(4)); %lamdadot2
ds = ds';
end

function con = bcfunc(so,sf)
% so is the initial states
% sf is the final states
con(1) = so(1); % x_1(0) = 0
con(2) = so(2); % x_2(0) = 0
con(3) = sf(2)-1; % x_2(tf) = 1
con(4) = sf(3)-10*sf(1); % Lambda_1(tf) = 10*x_1(tf)
end

function con = bcfunc2(so,sf)
% so is the initial states
% sf is the final states
%in the case of free time... the xo values are indeterminable, but the
%lambdao values are zero thus...
con(1) = so(3);% lambda1 = 0
con(2) = so(4);% lambda2 = 0
con(3) = sf(2)-1; %x_2(tf) = 1
con(4) = sf(3)-10*sf(1);% lambda1 = 10*x1(tf)
end