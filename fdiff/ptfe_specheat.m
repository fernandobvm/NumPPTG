function C = ptfe_specheat(T)
    %ptfe properties from Clark 1971
    C_r = 708.3;
    C_mc = 1220;
    T_m = 600;
    T_r = 298;
    C_h = 1537; 
    C_ma = 1476; %dont know if its equals to rho_m2
    T_h = 889;
    
    %cristaline phase
    if T < 600 %600K is maximum temperature that ptfe still cristaline
        C = C_r + (C_mc - C_r).*(T - T_r)/(T_m-T_r);
    %amorphous phase
    else
        C = C_h +(C_h - C_ma).*(T-T_m)/(T_h - T_m);
    end
end