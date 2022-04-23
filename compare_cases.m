clc; clear all; close all;

file1 = 'ppt_parallel_slug.mat';
file2 = 'ppt_parallel_slug.mat';
%file2 = 'ppt_parallel_modslug.mat';

load(file1);
v1 = voltage;
t1 = time;
i_bit1 = i_bit;
Isp1 = isp;
v_e1 = v_e;
n_t1 = n_t;
x_s1 = x_s;
xdot_s1 = xdot_s;
Ec1 = Ec;
Eb1 = Eb;
Ekinetic1 = Ekinetic;
Eohm1 = Eohm;

load(file2);
v2 = voltage;
t2 = time;
i_bit2 = i_bit;
Isp2 = isp;
v_e2 = v_e;
n_t2 = n_t;
x_s2 = x_s;
xdot_s2 = xdot_s;
Ec2 = Ec;
Eb2 = Eb;
Ekinetic2 = Ekinetic;
Eohm2 = Eohm;

% %comparativo tensão
% figure
% plot(t1,v1,'-o')
% hold on
% plot(t2,v2,'-o')
% title('Comparativo Tensão')
% hold off
% 
% %delta de tensão
% if length(t1) > length(t2)
%     v2_redimen = interp1(t2,v2,t1);
%     figure
%     plot(t1,v2_redimen-v1)
%     title('Delta da Tensão')
% else
%     v1_redimen = interp1(t1,v1,t2);
%     figure
%     plot(t2,v2-v1_redimen)
%     title('Delta da Tensão')
% end
% 
% %comparativo posição
% figure
% plot(t1,x_s1)
% hold on
% plot(t2,x_s2)
% title('Comparativo Posição')
% hold off
% 
% if length(t1) > length(t2)
%     x_s2_redimen = interp1(t2,x_s2,t1);
%     figure
%     plot(t1,x_s2_redimen-x_s1)
%     title('Delta da Posição')
% else
%     x_s1_redimen = interp1(t1,x_s1,t2);
%     figure
%     plot(t2,x_s2-x_s1_redimen)
%     title('Delta da Posição')
% end
% 
% %comparativo velocidade
% figure
% plot(t1,xdot_s1)
% hold on
% plot(t2,xdot_s2)
% title('Comparativo Velocidade')
% hold off
% 
% if length(t1) > length(t2)
%     xdot_s2_redimen = interp1(t2,xdot_s2,t1);
%     figure
%     plot(t1,xdot_s2_redimen-xdot_s1)
%     title('Delta da Velocidade')
% else
%     xdot_s1_redimen = interp1(t1,xdot_s1,t2);
%     figure
%     plot(t2,xdot_s2-xdot_s1_redimen)
%     title('Delta da Velocidade')
% end
% 
% i_bit = 100*abs(i_bit2/i_bit1 - 1)
% Isp = 100*abs(Isp2/Isp1 - 1)
% v_e = 100*abs(v_e2/v_e1 - 1)
% n_t = 100*abs(n_t2/n_t1 - 1)
% 


vondra = readtable('voltage.csv');
vondra.Properties.VariableNames = {'time','voltage'};
cap_energy = readtable('capenergy.csv');
cap_energy.Properties.VariableNames = {'time','cap_energy'};
ohm_energy = readtable('ohmenergy.csv');
ohm_energy.Properties.VariableNames = {'time','ohm_energy'};
ind_energy = readtable('indenergy.csv');
ind_energy.Properties.VariableNames = {'time','ind_energy'};
kin_energy = readtable('kinenergy.csv');
kin_energy.Properties.VariableNames = {'time','kin_energy'};
voltage2 = interp1(vondra.time, vondra.voltage,t1);
plot(t1,voltage2)
hold on
plot(t1,v1)
yticks([-500 0 500 1000 1500])
figure
plot(t1,v1-voltage2)
figure
plot(t1,(v1-voltage2)./voltage2)



figure
plot(vondra.time,vondra.voltage)
hold on
plot(t1,v1)
title('Comparativo Tensão no Capacitor')
xlabel('Tempo (s)')
ylabel('Tensão no capacitor (V)')
legend('Experimental','Numérico')
hold off

figure
plot(cap_energy.time,cap_energy.cap_energy)
hold on
plot(t1,Ec1)
title('Comparativo Energia no Capacitor')
xlabel('Tempo (s)')
ylabel('Energia no capacitor (J)')
legend('Experimental','Numérico')
hold off

figure
plot(ohm_energy.time,ohm_energy.ohm_energy)
hold on
plot(t1,Eohm1)
title('Comparativo Energia Dissipada')
xlabel('Tempo (s)')
ylabel('Energia dissipada (J)')
legend('Experimental','Numérico')
hold off

figure
plot(ind_energy.time,ind_energy.ind_energy)
hold on
plot(t1,Eb1)
title('Comparativo Energia no Indutor')
xlabel('Tempo (s)')
ylabel('Energia no indutor (J)')
legend('Experimental','Numérico')
hold off

figure
plot(kin_energy.time,kin_energy.kin_energy)
hold on
plot(t1,Ekinetic1)
title('Comparativo Energia Cinética')
xlabel('Tempo (s)')
ylabel('Energia cinética (J)')
legend('Experimental','Numérico')
hold off

