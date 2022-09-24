clc;clear all

n = 2e22;%n = 35*10^24;
k = 1.380648*10^-23;
m0 =1.66*10^-23;
mptfe = m0*(12.011*2 + 4*19.00)/6;
M = 100.0220/1000;
NA = 6.02e23;
m = mptfe;
Ts = 1000; %Ts = 700;
h = 0.03; w = 0.01; 
% h = 0.0254; w = h;
r_in = 0.75e-3; r_out = 3e-3;
A = h*w; 
% A = pi*(r_out^2-r_in^2);

atomic_emission_rate = sqrt(8*k*Ts/(pi*m))*n/4;

atomic_emission = atomic_emission_rate*A;

mol_emission = atomic_emission/NA;

mass_emission = mol_emission*M;

t_final = 2.5e-6; 
% t_final = 6.25e-6;
total_mass = t_final*mass_emission;
mass = t_final*atomic_emission*m

Pamb = 100 %10^5mbar
Pc = 1.847*10^15*Pamb;
Tc = 20815;
n0 = Pc*exp(-Tc/Ts)/(k*Ts)



Pv= n*k*Ts;
tal = Pv*sqrt(m/(2*pi*k*Ts));
mass2 = tal*A*t_final


