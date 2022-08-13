function [k, rho, C] = ptfe_properties(T)
    k_r = 0.2485;
    k_mc = 0.3619;
    T_m = 600;
    T_r = 298;
    k_h = 0.2472; 
    k_ma = 0.2453; %dont know if its equals to rho_m2
    T_h = 889;
    rho_r = 1933;
    rho_mc = 2174;
    rho_h = 1086; 
    rho_ma = 1740; %dont know if its equals to rho_m2
    C_r = 708.3;
    C_mc = 1220;
    C_h = 1537; 
    C_ma = 1476; %dont know if its equals to rho_m2
    
    %cristaline phase
    if T < 600 %600K is maximum temperature that ptfe still cristaline
        k = k_r + (k_mc - k_r).*(T - T_r)/(T_m-T_r);
        rho = rho_r + (rho_mc - rho_r).*(T - T_r)/(T_m-T_r);
        C = C_r + (C_mc - C_r).*(T - T_r)/(T_m-T_r);
    %amorphous phase
    else
        k = k_h +(k_h - k_ma).*(T-T_m)/(T_h - T_m);
        rho = rho_h +(rho_h - rho_ma).*(T-T_m)/(T_h - T_m);
        C = C_h +(C_h - C_ma).*(T-T_m)/(T_h - T_m);
    end
%     k = 0.2477 + 3.1717*10^-4*T;
%     rho = 1914;
%     C = 707.9;
end