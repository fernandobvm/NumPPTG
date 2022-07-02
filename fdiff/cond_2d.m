m = 40;
n = 50;

L = 1;
h = 0.2;

dx = L/(m-1);
dy = h/(n-1);

A = (2/dx^2 + 2/dy^2);
B = 1/dx^2;
C = 1/dy^2;


T = eye(m*n,m*n);
% T = zeros(m*n,m*n);
cc = zeros(m*n,1);

for i=2:n-1
    for j=2:m-1
%         T((n-1)*i,j+m) = 1;
        T(m*(i-1)+j,m*(i-1)+j) = 1;
        T(m*(i-1)+j,m*(i-1+1)+j) = -B/A;
        T(m*(i-1)+j,m*(i-1-1)+j) = -B/A;
        T(m*(i-1)+j,m*(i-1)+(j-1)) = -C/A;
        T(m*(i-1)+j,m*(i-1)+(j+1)) = -C/A;
    end
end

cc(1:m) = 1000; %superficie superior
cc(n*m-m+1:end) = 250;
x = T\cc;

fileID = fopen('temp_field.txt','w');
fprintf(fileID,'Variables="x","y","T"\n');
%fprintf(fileID,'Zone F=POINT,I=%d,J=%d\n\n',m,n);
result = zeros(m*n,3);
for i=1:n
    for j=1:m
        fprintf(fileID,'%f %f %f\n',dx*(j-1), h-dy*(i-1), x(m*(i-1)+j));
    end
end
k = 1;
for i=1:n
    for j=1:m
        result(k,:) = [dx*(j-1), h-dy*(i-1), x(m*(i-1)+j)];
        k = k+1;
    end
end
x = result(:,1);
y = result(:,2);
T = result(:,3);
z = zeros(length(x),1);
% vtkwrite('teste.vtk','structured_points', 'HeatTransfer', result)
vtkwrite('teste2.vtk','structured_grid',x,y,z,'scalars','Temperature',T)

