L = 1;
m = 15;
x = zeros(1,m);
x(end) = L;
for i=1:m-2
    x(1+i) = L/(2^i) + x(i);
end


m = 20;

x = zeros(1,m);
x(end) = L;

x1 = linspace(L/2,L,3*m/4+1);
x2 = linspace(0,L/2,m/4)
x = [x2(1:end-1), x1];