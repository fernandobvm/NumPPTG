clear all;
close all;
clc;

m = 80;

L = 1;
dx = L/(m-1);
% x = 0:dx:L;
x1 = linspace(3*L/5,L,1*m/8+1);
x2 = linspace(0,2*L/5,7*m/8);
x = [x2(1:end-1), x1];

dx = zeros(1,m-1);
for i=2:m
    dx(i-1) = x(i) - x(i-1);
end


% alfa = 0.124e-6;
% alfa = 9.7e-5;
alfa = 1e-4;

tf = 25;
steps = 100;
% dt = tf/steps;
dt = (min(dx))^2/(2*alfa);  %criterio de estabildiade - máximo dt
t = 0:dt/2:tf;


k = 0.3;
rho = 2.3e3;
% d = alfa*dt/dx^2;
% Fo = d;


T0 = 300;
q = 40000;
Twall = 0;
T = ones(m,length(t))*T0;
% T(1,:) = 0;
% T(m,:) = 0;

figure
for i= 2:length(t)
    for j=2:m-1
        d = alfa*(t(i)-t(i-1))/(x(j)-x(j-1))^2;
        T(j,i) = T(j,i-1) + d*(T(j-1,i-1) + T(j+1,i-1) - 2*T(j,i-1));
    end
%     T(1,i) = T(1,i-1) + d*(T(2,i-1) - 2*T(1,i-1));
    d = alfa*(t(i)-t(i-1))/(x(m)-x(m-1))^2;
    T(1,i) = Twall;
%     T(m,i) = T(m,i-1) + d*(T(m-1,i-1) - 2*T(m,i-1));
    plot(T(:,i))
    title(['t = ', num2str(t(i)), ' s'])
    pause(0.05)
end
T'

T_a = zeros(size(T));
% x = L - x;
for i=2:length(t)
        T_a(:,i) = (T0-Twall).*erf(x./(2*sqrt(alfa.*t(i)))) + Twall;
end
figure
plot(x,T_a(:,10))
hold on
scatter(x,T(:,10),'x')
xlabel('Comprimento do propelente (m)')
ylabel('Temperatura')
title('Perfil de Temperatura')
legend('Analitic', 'Numerical')

figure
p1 = plot(x,T_a(:,10),'b',x,T(:,10),'x')
hold on
p2 = plot(x,T_a(:,50),'r',x,T(:,50),'x')
p3 = plot(x,T_a(:,end),'r',x,T(:,end),'x')
xlabel('Comprimento do propelente (m)')
ylabel('Temperatura')
xlim([0 L/5])
title('Perfil de Temperatura')
legend([p1(1);p2(1);p3(1)],['t = ', num2str(t(15)), ' s'], ['t = ', num2str(t(50)), ' s'], ['t = ', num2str(t(end)), ' s'],'Location','southeast')



lT_a = T_a(:,10);
lT = T(:,10);
% save('data/m80_7-8_dt1-10_t=10.mat','lT_a','lT')