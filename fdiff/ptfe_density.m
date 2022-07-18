function rho = ptfe_density(T)
    %ptfe properties from Clark 1971
    rho_r = 1933;
    rho_mc = 2174;
    T_m = 600;
    T_r = 298;
    rho_h = 1086; 
    rho_ma = 1740; %dont know if its equals to rho_m2
    T_h = 889;
    
    %cristaline phase
    if T < 600 %600K is maximum temperature that ptfe still cristaline
        rho = rho_r + (rho_mc - rho_r).*(T - T_r)/(T_m-T_r);
    %amorphous phase
    else
        rho = rho_h +(rho_h - rho_ma).*(T-T_m)/(T_h - T_m);
    end
end