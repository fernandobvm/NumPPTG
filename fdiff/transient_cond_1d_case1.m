clear all;
close all;
clc;

tic;
m = 1000;

L = 1;
dx = L/(m-1);
% x = 0:dx:L;
x1 = linspace(4*L/5,L,7*m/8+1);
x2 = linspace(0,4*L/5,m/8);
x = [x2(1:end-1), x1];

dx = zeros(1,m-1);
for i=2:m
    dx(i-1) = x(i) - x(i-1);
end


% alfa = 0.124e-6;
% alfa = 9.7e-5;
alfa = 1e-4;

tf = 5;
steps = 100;
% dt = tf/steps;
dt = (min(dx))^2/(2*alfa);  %criterio de estabildiade - m�ximo dt
t = 0:dt/20:tf;


k = 0.3;
rho = 2.3e3;
% d = alfa*dt/dx^2;
% Fo = d;


T0 = 300;
q = 40000;

T = ones(m,length(t))*T0;
% T(1,:) = 0;
% T(m,:) = 0;

% figure
for i= 2:length(t)
    for j=2:m-1
        d = alfa*(t(i)-t(i-1))/(x(j)-x(j-1))^2;
        T(j,i) = T(j,i-1) + d*(T(j-1,i-1) + T(j+1,i-1) - 2*T(j,i-1));
    end
%     T(1,i) = T(1,i-1) + d*(T(2,i-1) - 2*T(1,i-1));
    d = alfa*(t(i)-t(i-1))/(x(m)-x(m-1))^2;
    T(m,i) = T(m,i-1) + d*(T(m-1,i-1) - 2*T(m,i-1) + T(m-1,i-1) + 2*(x(m)-x(m-1))*q/k);
%     plot(T(:,i))
%     title(['t = ', num2str(t(i)), ' s'])
%     pause(0.05)
end
% T'

T_a = zeros(size(T));
x = L - x;
for i=2:length(t)
        T_a(:,i) = (q./k).*sqrt(4*alfa.*t(i)/pi).*(exp(-x.^2./(4*alfa.*t(i))) - x.*sqrt(pi./(4*alfa.*t(i))).*erfc(x./(sqrt(4*alfa.*t(i))))) + T0;
end
figure
plot(L - x,T_a(:,end))
hold on
scatter(L - x,T(:,end),'x')
xlabel('Comprimento do propelente (m)')
ylabel('Temperatura')
title('Perfil de Temperatura')
legend('Analitic', 'Numerical')

toc
lT_a = T_a(:,10);
lT = T(:,10);
% save('data/m80_7-8_dt1-10_t=10.mat','lT_a','lT')