
tf = 5;
dt = (min(dx))^2/(2*alfa);  %criterio de estabildiade - máximo dt
t = 0:dt/2:tf;

T0 = 1;
m = 800;
T = ones(m, length(t))*T0;

for i= 2:length(t)
    for j=2:m-1
        d = alfa*(t(i)-t(i-1))/(x(j)-x(j-1))^2;
        T(j,i) = T(j,i-1) + d*(T(j-1,i-1) + T(j+1,i-1) - 2*T(j,i-1));
    end
    d = alfa*(t(i)-t(i-1))/(x(m)-x(m-1))^2;
    T(m,i) = T(m,i-1) + d*(T(m-1,i-1) - 2*T(m,i-1) + T(m-1,i-1) + 2*(x(m)-x(m-1))*q/k);
end
