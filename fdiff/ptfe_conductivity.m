function k = ptfe_conductivity(T)
    %ptfe properties from Clark 1971
    k_r = 0.2485;
    k_mc = 0.3619;
    T_m = 600;
    T_r = 298;
    k_h = 0.2472; 
    k_ma = 0.2453; %dont know if its equals to rho_m2
    T_h = 889;
    
    %cristaline phase
    if T < 600 %600K is maximum temperature that ptfe still cristaline
        k = k_r + (k_mc - k_r).*(T - T_r)/(T_m-T_r);
    %amorphous phase
    else
        k = k_h +(k_h - k_ma).*(T-T_m)/(T_h - T_m);
    end
%     k = 0.2477 + T*3.1717*10^-4;
end