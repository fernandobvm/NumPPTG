clc;
clear all;
close all;

load('data/m80_3-4_dt1-2_t=10')
T_a1 = lT_a;
T1 = lT;

load('data/m80_7-8_dt1-2_t=10')
T_a2 = lT_a;
T2 = lT;

load('data/m160_3-4_dt1-2_t=10')
T_a3 = lT_a;
T3 = lT;

load('data/m160_7-8_dt1-2_t=10')
T_a4 = lT_a;
T4 = lT;

clear lT_a lT;

plot(T_a1(50:end))
hold on
plot(T1(50:end))
plot(T2(50:end))
plot(T3(130:end))
plot(T4(130:end))
plot(T_a4(130:end))
legend('Analitic', 'Numerical1','Numerical2','Numerical3','Numerical4','Analitic4')

load('data/m80_3-4_dt1-10_t=10')
T_a1 = lT_a;
T1 = lT;

load('data/m80_7-8_dt1-10_t=10')
T_a2 = lT_a;
T2 = lT;

load('data/m160_3-4_dt1-10_t=10')
T_a3 = lT_a;
T3 = lT;

load('data/m160_7-8_dt1-10_t=10')
T_a4 = lT_a;
T4 = lT;

clear lT_a lT;
figure
plot(T_a1(50:end))
hold on
plot(T1(50:end))
plot(T2(50:end))
plot(T3(130:end))
plot(T4(130:end))
plot(T_a4(130:end))
legend('Analitic', 'Numerical1','Numerical2','Numerical3','Numerical4','Analitic4')