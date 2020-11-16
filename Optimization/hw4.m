clear
close all
clc

x_0 = 0.75;

t_f = -log(1-x_0);

t = linspace(0,t_f);

x = 1-exp(t-t_f);

figure(1)

plot(t, x, '--b');
xlabel('time, (s)')
ylabel('x, (m)')

syms t
nu = 2*exp(t-2);
x_o = .75
eqn = 1+(x_o-1)*exp(t) == -nu/4*(exp(t-2)-exp(2-t));

t1 = solve(eqn,t);

nu = double(subs(nu,t,t1(2)))

t1 = double(t1(2))

clear

J = 1.6217 - (1.37^2)/8*(1-exp(2*(2-1.6217)))


%% Numerical solutions
tf = 2;
tflip = linspace(0,tf,tf*100+1);

options = odeset('Events',@satEvent);
[T,S] = ode45(@(t,s)eom(t,s),tflip(end:-1:1),[0,1.37],options);
T = flip(T);
S = flip(S);

%[T1,S1] = ode45(@(t,s)eom(t,s),linspace()),[0,1.37],options);

figure(2)
plot(T,S(:,1),'--r')
hold on
% solinit.x = T;
% solinit.y = S';
% sol = bvp4c(@(tx,sx)ELeom(tx,sx),@bcfunc,solinit);
% plot(sol.x,sol.y(1,:),'-k');
% 
t = linspace(0,tf,tf*100+1);
xsat = (.75-1).*exp(t(1:163))+1;
xa = 1.37/4.*(exp(t(164:201)-2)-exp(2-t(164:201)));
xg = [xsat xa];

plot(t,xg,'ob');

hold off
% Euler-Lagrange Equations:
function ds = eom(~,s)
% s(1) = x, s(2) = Lambda
u = s(2)/2; % optimal control law
ds(1) = -(s(1)-u); %xdot = x-u
ds(2) = s(2); %lambda_dot = -lambda
ds = ds'; 
end

% Euler-Lagrange Equations:
function ds = ELeom(~,sx)
% sx(1) = x, sx(2) = Lambda
if sx(2)/2 > 1
    u = 1;
elseif sx(2)/2 < -1
    u = -1;
else
    u = sx(2)/2; % optimal control law
end
ds(1) = sx(1)-u; %xdot = x-u
ds(2) = -sx(2); %lambda_dot = -lambda
ds = ds';
end

% Boundary Condition Function
function con = bcfunc(so,sf)
% so is the initial states
% sf is the final states
%con(1) = so(1)-.75; % x_1(0) = 0
con(2) = sf(2)-1.37; % x_2(0) = 0
con(1) = sf(1);
end

function [position,isterminal,direction] = satEvent(t,s)
position = s(2)/2-1; %when u becomes saturated
direction = 0;
isterminal = 1; %halting integration
end