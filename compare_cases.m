%file1 = 'ppt_parallel_slug.mat';
file1 = 'ppt_parallel_slug1.mat';
file2 = 'ppt_parallel_modslug.mat';

load(file1);
v1 = voltage;
t1 = time;
i_bit1 = i_bit;
Isp1 = isp;
v_e1 = v_e;
n_t1 = n_t;
x_s1 = x_s;
xdot_s1 = xdot_s;

load(file2);
v2 = voltage;
t2 = time;
i_bit2 = i_bit;
Isp2 = isp;
v_e2 = v_e;
n_t2 = n_t;
x_s2 = x_s;
xdot_s2 = xdot_s;

%comparativo tensão
figure
plot(t1,v1)
hold on
plot(t2,v2)
title('Comparativo Tensão')
hold off

%delta de tensão
figure
plot(t1,v2-v1)
title('Delta da Tensão')

%comparativo posição
figure
plot(t1,x_s1)
hold on
plot(t2,x_s2)
title('Comparativo Posição')
hold off

figure
plot(t1,x_s2-x_s1)
title('Delta da Posição')

%comparativo velocidade
figure
plot(t1,xdot_s1)
hold on
plot(t2,xdot_s2)
title('Comparativo Velocidade')
hold off

figure
plot(t1,xdot_s2-xdot_s1)
title('Delta da Velocidade')


i_bit = 100*abs(i_bit2/i_bit1 - 1)
Isp = 100*abs(Isp2/Isp1 - 1)
v_e = 100*abs(v_e2/v_e1 - 1)
n_t = 100*abs(n_t2/n_t1 - 1)
