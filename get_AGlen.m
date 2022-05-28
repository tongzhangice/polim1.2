function [AGlen_s] = get_AGlen(T, omega, CTS, para)
% Date: 2017-12-17
% Author: Wang Yuzhe
% Calculate the flow rate factor A.

% Input:
% T: temperature
% para: parameters (struct)

% Output:
% AGlen_s: flow rate factor A on staggered grid.

global N M H zeta

SPY = para.SPY;
rhoi = para.rhoi;
g = para.g;
R = para.R;
betaCC = para.betaCC;
type_Arrhenius = para.type_Arrhenius;
is_duval = para.is_duval;

AGlen = zeros(N,M);
if strcmpi(type_Arrhenius,'greve')
    AGlen0_cold = 3.985e-13;
    AGlen0_warm = 1.916e3;
    
    Q_cold = 60e3;
    Q_warm = 139e3;
    for i = 1:M
        T_Corr = T(:,i) + betaCC*rhoi*g*(1-zeta).*H(i);
        idx_le_n10 = (T_Corr <= 263.15);
        idx_ge_n10 = (T_Corr > 263.15);
        
        AGlen(:,i) = AGlen0_cold*exp(-Q_cold./(R*T_Corr)).*idx_le_n10 +...
            AGlen0_warm*exp(-Q_warm./(R*T_Corr)).*idx_ge_n10; % [Pa-3 s-1]
    end
    AGlen = AGlen*SPY; % [Pa-3 a-1]
    
elseif strcmpi(type_Arrhenius,'cuffey')
    AGlen0 = 3.5e-25; % constant prefactor at -10 celsius [Pa-3 s-1], Cuffey&Paterson (2010), p74, eq3.36
    Q_cold = 60e3; % activation energy for creep [J mol-1]
    Q_warm = 115e3; % activation energy for creep [J mol-1]
    
    for i = 1:M
        delta_Tm = betaCC*rhoi*g*H(i)*(1-zeta);
        Tref_Corr = 263.15 + delta_Tm;
        T_Corr = T(:,i) + delta_Tm;
        
        % cold layer
        idx_cold = (CTS(i)+1):N;
        idx_le_n10 = (T(idx_cold,i)<=263.15);
        idx_ge_n10 = (T(idx_cold,i)>263.15);
        T_Corr_cold = T_Corr(idx_cold);
        Tref_Corr_cold = Tref_Corr(idx_cold);
        AGlen(idx_cold,i) = AGlen0*exp(-Q_cold./R.*(1./T_Corr_cold-1./Tref_Corr_cold)).*idx_le_n10 +...
            AGlen0*exp(-Q_warm./R.*(1./T_Corr_cold-1./Tref_Corr_cold)).*idx_ge_n10; % [Pa-3 s-1]
        
        % temperate layer
        if CTS(i)~=0
            idx_temp = 1:CTS(i);
            idx_O_le = (omega(idx_temp,i)<=para.moi_UB); % moisture less than upper_bound
            idx_O_ge = (omega(idx_temp,i)>para.moi_UB); % moisture larger than upper_bound
            omega(idx_temp,i) = omega(idx_temp,i).*idx_O_le + para.moi_UB*idx_O_ge;
            
            % Duval-Lliboutry relation for temperate ice layer
            AGlen(idx_temp,i) = 24e-25*(1 + is_duval*181.25*omega(idx_temp,i));
        end
    end
    AGlen = AGlen*SPY; % [Pa-3 a-1]
    
else
    disp('type_Arrhenius should be "greve", "cuffey"')
end

AGlen_s = main2staggerX(AGlen);
end