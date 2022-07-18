rho_r = 1933;
rho_mc = 2174;
T_m = 600;
T_r = 298;

T_c = [0 298 300 311 317 317.5 322 345 450 554 582 599];
rho_c = rho_r + (rho_mc - rho_r).*(T_c - T_r)/(T_m-T_r);

rho_h = 1086; 
rho_ma = 1740; %dont know if its equals to rho_m2
T_h = 889;

T_a = [601 609 623 689 710 715 769 781 800 900 1000];
rho_a = rho_h +(rho_h - rho_ma).*(T_a-T_m)/(T_h - T_m);

T = [T_c T_a];
rho = [rho_c rho_a];
plot(T,rho)