% clc;
% clear all;
% 
% L = 100;
% N = 100;
% dx = L/(N-1);
% 
% T = zeros(N,1);
% 
% Tb = 100;
% 
% 
% k = 100; %iterations
% 
% for j=1:1:k
%     T(1,1) = Tb;
%     for i = 2:1:N-1
%         T(i,1) = (T(i+1,1)+T(i-1,1))/2;
%     end
%     T(N,1) = T(N-1,1);
% end
% 
% plot(T)

L = 0.1;
n = 10;
T0 = 0;
Te = 40;
Td = 20;

dx = L/n;
alpha = 0.0001;

tf = 60;
dt = 0.1;

x = dx/2:dx:L-dx/2;

t = 0:dt:tf;

T = ones(n,length(t))*T0;
dTdt = zeros(n,1);

% for j=2:length(t)
%     for i=2:n-1
%         dTdt(i) = alpha*(-(T(i,j-1)-T(i-1,j-1))/dx^2 + (T(i+1,j-1)-T(i,j-1))/dx^2);
%     end
%     dTdt(1) = alpha*(-(T(1,j-1)-Te)/dx^2 + (T(2,j-1)-T(1,j-1))/dx^2); 
%     dTdt(n) = alpha*(-(T(n,j-1)-T(n-1,j-1))/dx^2 + (Td-T(n,j-1))/dx^2);
%     T(:,j) = T(:,j) + dTdt*dt;
% end

for j=1:length(t)
    for i=2:n-1
        dTdt(i) = alpha*(-(T(i)-T(i-1))/dx^2 + (T(i+1)-T(i))/dx^2);
    end
    dTdt(1) = alpha*(-(T(1)-Te)/dx^2 + (T(2)-T(1))/dx^2); 
    dTdt(n) = alpha*(-(T(n)-T(n-1))/dx^2 + (Td-T(n))/dx^2);
    T = T + dTdt*dt;
    
    figure(1)
    plot(x,T)
    axis([0 L 0 50])
%     pause(0.1)
end
