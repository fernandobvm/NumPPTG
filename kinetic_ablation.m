T0 = 900; %K
T2 = 23209.05; %2ev
n2 = 1e21;
Pc = 1.84*10^15; %N/m^2
Tc = 20815; %K
k = 1.380648*10^-23;
% M = 100.16/1000;
% NA = 6.0221408e+23;
% m = M/NA;
m0 =1.66*10^-23;
mptfe = m0*(12.011*2 + 4*19.00)/6;
m = mptfe;

syms n0 n1 T1 V1 beta

eq1 = n0/(2*sqrt(pi*m/2*k*T0)) == n1*V1 + beta*n1*(exp(-(V1/sqrt(2*k*T1/m))^2)-(V1/sqrt(2*k*T1/m))*sqrt(pi)*erfc((V1/sqrt(2*k*T1/m))))/(2*sqrt(pi*m/2*k*T1)) ;
eq2 = n0*2*k*T0/(4*m) == n1*((1+2*(V1/sqrt(2*k*T1/m))^2) - beta*((0.5 + (V1/sqrt(2*k*T1/m))^2)*erfc((V1/sqrt(2*k*T1/m))) - (V1/sqrt(2*k*T1/m))*exp(-(V1/sqrt(2*k*T1/m))^2)/sqrt(pi)))*2*k*T1/(2*m);
eq3 = n0/(pi*m/(2*k*T0))^1.5 == (V1/sqrt(2*k*T1/m))*n1*(((V1/sqrt(2*k*T1/m))^2 + 2.5) - beta*(((V1/sqrt(2*k*T1/m))^2+2.5)*erfc((V1/sqrt(2*k*T1/m))) - (2 + (V1/sqrt(2*k*T1/m))^2)*exp(-(V1/sqrt(2*k*T1/m))^2)/((V1/sqrt(2*k*T1/m))*sqrt(pi)))/2)/(pi*(m/(2*k*T1))^1.5);
eq4 = (V1/sqrt(2*k*T1/m))^2 == (n2*T2/(2*T1) - n1/2)/(n1*(1 - n1/n2));
eq5 = Pc*exp(-Tc/T0) == n0*k*T0;
eqns = [eq1, eq2, eq3, eq4, eq5];

%S = solve(eqns);

%-----------------------------------------------------------------------
%ALL EQUATIONS
% eq1 = n0/(2*sqrt(pi*d0)) == n1*V1 + ((beta*n1)/(2*sqrt(pi*d1)))*(exp(-alfa^2) - alfa*sqrt(pi)*erfc(alfa)); %OKAY
% eq2 = n0/(4*d0) == (n1/(2*d1))*((1 + 2*alfa^2) - beta*((0.5 + alfa^2)*erfc(alfa) - alfa*exp(-alfa^2)/sqrt(pi))); %OKAY
% eq3 = n0/(pi*d0)^1.5 == (n1/(pi*d1^1.5))*(alfa*(alfa^2 + 2.5) - 0.5*beta*((2.5 + alfa^2)*alfa*erfc(alfa) - (2 + alfa^2)*exp(-alfa^2)/sqrt(pi))); %OKAY
% eq4 = d0 == m/(2*k*T0); %OKAY
% eq5 = d1 == m/(2*k*T1);%OKAY
% eq6 = V1^2/(2*k*T1/m) == ((T2*n2/(2*T1)) - (n1/2))/(n1*(1 - n1/n2));%OKAY
% eq7 = Pc*exp(-Tc/T0) == n0*k*T0;%OKAY

syms n0 n1 T1 V1 beta
eq1 = n0/(2*sqrt(pi*(m/(2*k*T0)))) == n1*V1 + ((beta*n1)/(2*sqrt(pi*(m/(2*k*T1)))))*(exp(-(V1^2/(2*k*T1/m))^2) - (V1^2/(2*k*T1/m))*sqrt(pi)*erfc((V1^2/(2*k*T1/m))));
eq2 = n0/(4*(m/(2*k*T0))) == (n1/(2*(m/(2*k*T1))))*((1 + 2*(V1^2/(2*k*T1/m))^2) - beta*((0.5 + (V1^2/(2*k*T1/m))^2)*erfc((V1^2/(2*k*T1/m))) - (V1^2/(2*k*T1/m))*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq3 = n0/(pi*(m/(2*k*T0)))^1.5 == (n1/(pi*(m/(2*k*T1))^1.5))*((V1^2/(2*k*T1/m))*((V1^2/(2*k*T1/m))^2 + 2.5) - 0.5*beta*((2.5 + (V1^2/(2*k*T1/m))^2)*(V1^2/(2*k*T1/m))*erfc((V1^2/(2*k*T1/m))) - (2 + (V1^2/(2*k*T1/m))^2)*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq4 = V1^2/(2*k*T1/m) == ((T2*n2/(2*T1)) - (n1/2))/(n1*(1 - n1/n2));
eq5 = Pc*exp(-Tc/T0) == n0*k*T0;
eqn = [eq1, eq2, eq3, eq4, eq5];
solution = solve(eqn)
%-----------------------------------------------------------------------




syms a b c d T1 V1 beta
eq1 = (a*10^b)/(2*sqrt(pi*(m/(2*k*T0)))) == (c*10^d)*V1 + ((beta*(c*10^d))/(2*sqrt(pi*(m/(2*k*T1)))))*(exp(-(V1^2/(2*k*T1/m))^2) - (V1^2/(2*k*T1/m))*sqrt(pi)*erfc((V1^2/(2*k*T1/m))));
eq2 = (a*10^b)/(4*(m/(2*k*T0))) == ((c*10^d)/(2*(m/(2*k*T1))))*((1 + 2*(V1^2/(2*k*T1/m))^2) - beta*((0.5 + (V1^2/(2*k*T1/m))^2)*erfc((V1^2/(2*k*T1/m))) - (V1^2/(2*k*T1/m))*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq3 = (a*10^b)/(pi*(m/(2*k*T0)))^1.5 == ((c*10^d)/(pi*(m/(2*k*T1))^1.5))*((V1^2/(2*k*T1/m))*((V1^2/(2*k*T1/m))^2 + 2.5) - 0.5*beta*((2.5 + (V1^2/(2*k*T1/m))^2)*(V1^2/(2*k*T1/m))*erfc((V1^2/(2*k*T1/m))) - (2 + (V1^2/(2*k*T1/m))^2)*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq4 = V1^2/(2*k*T1/m) == ((T2*n2/(2*T1)) - ((c*10^d)/2))/((c*10^d)*(1 - (c*10^d)/n2));
eq5 = Pc*exp(-Tc/T0) == (a*10^b)*k*T0;
eqn2 = [eq1, eq2, eq3, eq4, eq5];
solution = solve(eqn2)

%------------------------------------------

syms n1 T1 V1 beta
n0 = Pc*exp(-Tc/T0)/(k*T0);
eq1 = n0/(2*sqrt(pi*(m/(2*k*T0)))) == n1*V1 + ((beta*n1)/(2*sqrt(pi*(m/(2*k*T1)))))*(exp(-(V1^2/(2*k*T1/m))^2) - (V1^2/(2*k*T1/m))*sqrt(pi)*erfc((V1^2/(2*k*T1/m))));
eq2 = n0/(4*(m/(2*k*T0))) == (n1/(2*(m/(2*k*T1))))*((1 + 2*(V1^2/(2*k*T1/m))^2) - beta*((0.5 + (V1^2/(2*k*T1/m))^2)*erfc((V1^2/(2*k*T1/m))) - (V1^2/(2*k*T1/m))*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq3 = n0/(pi*(m/(2*k*T0)))^1.5 == (n1/(pi*(m/(2*k*T1))^1.5))*((V1^2/(2*k*T1/m))*((V1^2/(2*k*T1/m))^2 + 2.5) - 0.5*beta*((2.5 + (V1^2/(2*k*T1/m))^2)*(V1^2/(2*k*T1/m))*erfc((V1^2/(2*k*T1/m))) - (2 + (V1^2/(2*k*T1/m))^2)*exp(-(V1^2/(2*k*T1/m))^2)/sqrt(pi)));
eq4 = V1^2/(2*k*T1/m) == ((T2*n2/(2*T1)) - (n1/2))/(n1*(1 - n1/n2));
eq5 = Pc*exp(-Tc/T0) == n0*k*T0;
eqn = [eq1, eq2, eq3, eq4]%, eq5];
solution = solve(eqn)
