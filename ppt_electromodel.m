%% LES 6 PPT Data
clc;clear all;

Rc = 30e-3;
Re = 0; %não achei explicitamente declarado, teste zerado
Rpe =  0; %não achei explicitamente declarado, teste zerado
Rp = 0.074;
Req = Rc + Re + Rpe + Rp;
Lc = 0; 
%Lpe = %varia com o tempo
Le = 34e-9; %teste
L0 = 34e-9; %Lc+Le
C = 2e-6;
Te = 1.5;
n_e = 1e21;
t_pulse = 0.4e-6;

V0 = 1360;
h = 3e-2;
w = 1e-2;
l = 0.6e-2;
delta = 0; %thin current sheet
alpha = 1; %mass distribution model, alpha = 0 (uniform) alpha = 1 (slug model)
m_bit = 10e-6;
mu = 1.2566e-6;

tspan = [0 2.5e-6];
y0 = [0,0,0,0];
MaxStep = 1e-8;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep', MaxStep);

global iter
global teste
teste = -1*ones(1000000,1);

iter = 1;
%[t,y] = ode45(@(t,y) ppt(t,y,V0,C,Req,mu,h,w,m_bit,L0,delta,MaxStep), tspan, y0, opts);
[t,y] = ode45(@(t,y) ppt(t,y,V0,C,Req,mu,h,w,m_bit,Lc,Le,delta), tspan, y0, opts);

I = y(:,4);
V = V0 - y(:,2)/C;
Ec = C*V.^2/2;

figure
yyaxis left
plot(t,V)
yyaxis right
plot(t,I)

Eohm = zeros(length(t),1);
for i = 2:length(t)
    int = trapz(t(1:i),Req*y(1:i,4).*y(1:i,4));
    Eohm(i) = int;
end

figure
% Eohm = MaxStep*cumtrapz(Req*y(:,4).*y(:,4))/4;
%Eohm = trapz(t(1:5),Req*y(1:5,4).*y(1:5,4));
plot(t,Eohm)
hold on
plot(t,Ec)

function y = ppt1(t,y0,V0,C,Req,mu,h,w,m_bit,L0,delta,MaxStep)
    global iter
    global teste
    y = zeros(4,1);
    teste(iter) = y0(1);
    iter = iter + 1;
    y(1) = y0(2);
    integral = MaxStep*cumtrapz(teste(1:iter-1));
    y(2) = (V0 - integral(end)/C - y0(1)*Req - mu*h*y0(4)*y0(1)/w)/(L0+mu*h*y0(3)/w + mu*delta*h/(2*w));
    %y(2) = V0 - y0(1)/C - y0(1)*Req - mu*h*y0(4)*y0(1)/w;
    y(3) = y0(4);
    y(4) = mu*h*y(1)*y(1)/(2*m_bit*w);

end

function y = ppt(t,y0,V0,C,Req,mu,h,w,m_bit,Lc,Le,delta)
    Lpe = mu*h*y0(1)/w +mu*delta*h/(2*w);
    L = Lpe + Lc + Le;
    
    y = zeros(4,1);
    y(1) = y0(3);
    y(2) = y0(4);
    y(3) = Lpe*y0(4)*y0(4)/(2*m_bit);
    y(4) = (-y0(2)/C - mu*h*y0(3)*y0(4)/w - Req*y0(4) + V0)/(L);

end