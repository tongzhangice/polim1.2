function [ output_args ] = get_basalMeltRate(u_s, u, visc_s, T, para)
%GET_BASALMELTRATE Summary of this function goes here
%   Detailed explanation goes here

[tauxz, ~] = calc_tauxz(u_s, visc_s);

m_basal = (para.Qgeo +...
    para.kc*(T(2,:)-T(1,:))./(H*dzeta) +...
    tauxz.*u(1,:)/para.SPY) /(para.rhow*para.Lw)*para.SPY; % basal melt [m a-1]

m_basal(m_basal<0)=0;


end

