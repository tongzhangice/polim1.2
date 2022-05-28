% Date: 2017-7-12
% Author: Wang Yuzhe
% Calculate the glacier thermal regime using the enthalpy scheme.
% Use the enthalpy definition and the water drainage model proposed by Anderas Ascthk_wanden.

% <At_isTEMP>, <At_thk_TEMP>, and <At_thk_w> determine the basal boundary condition

function [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_greveDrain, qw_TEMP_diffu] = ...
    solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt, Esbc, Eini, para)
global M N dx dzeta zeta H dzetadx iTimeStep...
    At_E At_Kappa_vs At_omega...
    At_isTEMP At_thk_TEMP At_thk_w

SPY = para.SPY;
rhoi = para.rhoi;
rhow = para.rhow;
g = para.g;
kc = para.kc;
Cp = para.Cp;
Qgeo = para.Qgeo;
betaCC = para.betaCC;
Lw = para.Lw;
Tref = para.Tref;
Kc = para.Kc;
Kt = para.Kt;

is_auto_enth_BBC = para.is_auto_enth_BBC;
type_enth_BBC = para.type_enth_BBC;
has_Greve_drainage = para.has_Greve_drainage;

% convert [xx a-1] to [xx s-1]
dt = dt*SPY;
u = u/SPY; % [m s-1]
u_s = u_s/SPY; % [m s-1]
w = w/SPY;  % [m s-1]
w_vs = w_vs/SPY; % [m s-1]
strainHeat = strainHeat/SPY; % [Pa s-1];

% 'var_vs' means the secondary grid point in zeta coordinate
dzetadx_vs = (dzetadx(1:end-1,:) + dzetadx(2:end,:))/2; % <N-1 * M>
u_vs = (u(1:end-1,:) + u(2:end,:))/2; % <N-1 * M>

% coefficient for the vertical derivative of enthalpy (\partial E / \partial zeta)
coeff = u.*dzetadx + w./(ones(N,1)*H); % [s-1]; <N * M>
coeff_vs = u_vs.*dzetadx_vs + w_vs(1:N-1, :)./(ones(N-1,1)*H); % [s-1]; <N-1 * M>

% initialization
LT1 = zeros(N,1);
LT2 = zeros(N,1);
LT3 = zeros(N,1);
RT = zeros(N,1);

E = zeros(N,M);
T = zeros(N,M);
omega = zeros(N,M);

thk_TEMP = zeros(1,M);
m_basal = zeros(1,M);
Kappa_i = zeros(N,1);
Kappa_vs = zeros(N-1,M);
CTS = zeros(1,M);
qw_greveDrain = zeros(1,M);
drainageRate = zeros(N,1);
qw_TEMP_diffu = zeros(N-1,M);

if iTimeStep==1
    is_trans = 0;
    Enm1 = Eini;
    Kappa_vs_lst = Kc*ones(N-1,M);
    omega_lst = zeros(N,M);
else
    is_trans = 1;
    Enm1 = At_E(:,:,iTimeStep-1);
    Kappa_vs_lst = At_Kappa_vs(:,:,iTimeStep-1);
    omega_lst = At_omega(:,:,iTimeStep-1);
end

for i=1:M
    % Enthalpy at the pressure-melting point for the i-th column
    Tpmp_i = 273.15 - betaCC*rhoi*g*H(i)*(1-zeta); % <N * 1>
    Epmp_i = Cp*(Tpmp_i - Tref); % <N * 1>
    
    for j = 2:N-1
        if has_Greve_drainage
            drainageRate(j) = Greve_drain(omega_lst(j,i))/SPY; % [s-1]
        else
            drainageRate(j) = 0.0;
        end
        
        if M==1 % ice slab (no horizontal advection and diffusion)
            RT(j) = is_trans*Enm1(j,i)/(dt) + strainHeat(j,i)/rhoi -...
                rhow/rhoi*Lw*drainageRate(j);
        else
            if i==1
                RT(j) = is_trans*Enm1(j,i)/(dt) + strainHeat(j,i)/rhoi;
            else
                RT(j) = is_trans*Enm1(j,i)/(dt) + strainHeat(j,i)/rhoi...
                    + u_s(j,i)*Enm1(j,i-1)/dx -...
                    rhow/rhoi*Lw*drainageRate(j);
            end
        end
        
        if coeff(j,i)>0
            LT1(j) = -coeff_vs(j-1,i)/dzeta - Kappa_vs_lst(j-1,i)/(H(i)^2*dzeta^2);
            LT2(j) = is_trans/(dt) + u_s(j,i)/dx + coeff_vs(j-1,i)/dzeta +...
                (Kappa_vs_lst(j-1,i)+Kappa_vs_lst(j,i))/(H(i)^2*dzeta^2);
            LT3(j) = -Kappa_vs_lst(j,i)/(H(i)^2*dzeta^2);
        elseif coeff(j,i)<=0
            LT1(j) = -Kappa_vs_lst(j-1,i)/(H(i)^2*dzeta^2);
            LT2(j) = is_trans/(dt) + u_s(j,i)/dx - coeff_vs(j,i)/dzeta +...
                (Kappa_vs_lst(j-1,i)+Kappa_vs_lst(j,i))/(H(i)^2*dzeta^2);
            LT3(j) = coeff_vs(j,i)/dzeta - Kappa_vs_lst(j,i)/(H(i)^2*dzeta^2);
        end
    end
    
    % surface BC
    LT1(N) = 0;
    LT2(N) = 1;
    LT3(N) = 0;
    RT(N) = Esbc(i);
    
    % basal BC
    if is_auto_enth_BBC % decision chart (Aschwanden et al., 2012, Figure 5)
        if iTimeStep==1
            [LT, RT] = enth_BBC_ColdBaseDry(LT1, LT2, LT3, RT, H(i), para);
        else
            if At_isTEMP(iTimeStep-1,i)==0
                if At_thk_w(iTimeStep-1,i) > 0
                    [LT, RT] = enth_BBC_ColdBaseWet(LT1, LT2, LT3, RT, Epmp_i);
                else
                    [LT, RT] = enth_BBC_ColdBaseDry(LT1, LT2, LT3, RT, H(i), para);
                end
            else
                if At_thk_TEMP(iTimeStep-1,i) > 0
                    [LT, RT] = enth_BBC_TempLayer(LT1, LT2, LT3, RT);
                else
                    [LT, RT] = enth_BBC_TempBase(LT1, LT2, LT3, RT, Epmp_i);
                end
            end
        end
    else % Specify a type of basal boundary condition
        switch type_enth_BBC
            case 1
                [LT, RT] = enth_BBC_ColdBaseDry(LT1, LT2, LT3, RT, H(i), para);
            case 2
                [LT, RT] = enth_BBC_TempBase(LT1, LT2, LT3, RT, Epmp_i);
            case 3
                [LT, RT] = enth_BBC_TempLayer(LT1, LT2, LT3, RT);
            case 4
                [LT, RT] = enth_BBC_ColdBaseWet(LT1, LT2, LT3, RT, Epmp_i);
        end
    end
    
    % solution
    E(:,i) = LT\RT; % <N * 1>
    
    logic1 = E(:,i) >= Epmp_i; % 1: temperate; 0: cold. <N * 1>
    T(:,i) = logic1.*Tpmp_i + (1-logic1).*(E(:,i)/Cp + Tref);
    omega(:,i) = logic1.*(E(:,i)-Epmp_i)/Lw;    
    jcts = find(logic1, 1, 'last'); % empty: cold; jcts>=1: temperate
    
    % corrector step (Blatter & Greve, 2015, p201)
    if isempty(jcts)      %% cold
        CTS(i) = 0;
        thk_TEMP(i) = 0;
        m_basal(i) = 0;
        Kappa_vs(:,i) = Kappa_vs_lst(:,i); %Kc*ones(N-1,1);
    else                  %% polythermal
        CTS(i) = jcts;
        thk_TEMP(i) = (jcts-1)*H(i)*dzeta;
        %         Kappa_vs(:,i) = Kappa_i(1:end-1).^(0.5).*Kappa_i(2:end).^(0.5); % geomtetric mean
        Kappa_i(1:jcts-1) = Kt;
        Kappa_i(jcts:end) = Kc;
        Kappa_harmmean = harmmean([Kappa_i(1:end-1)'; Kappa_i(2:end)']);
        Kappa_vs(:,i) = Kappa_harmmean'; % harmonic mean

        % exterior boundary condition
        if is_auto_enth_BBC
            if jcts==1
                [LT, RT] = enth_BBC_TempBase(LT1, LT2, LT3, RT, Epmp_i);
            else
                [LT, RT] = enth_BBC_TempLayer(LT1, LT2, LT3, RT);
            end
        else
            switch type_enth_BBC
                case 1
                    [LT, RT] = enth_BBC_ColdBaseDry(LT1, LT2, LT3, RT, H(i), para);
                case 2
                    [LT, RT] = enth_BBC_TempBase(LT1, LT2, LT3, RT, Epmp_i);
                case 3
                    [LT, RT] = enth_BBC_TempLayer(LT1, LT2, LT3, RT);
                case 4
                    [LT, RT] = enth_BBC_ColdBaseWet(LT1, LT2, LT3, RT, Epmp_i);
            end
        end
        
        % update enthalpy
        E(:,i) = LT\RT; % <N * 1>

        % inverted temperature and water content
        logic1 = E(:,i)>=Epmp_i; % 1: temperate; 0: cold. <N * 1>
        T(:,i) = logic1.*Tpmp_i + (1-logic1).*(E(:,i)/Cp + Tref);
        omega(:,i) = logic1.*(E(:,i)-Epmp_i)/Lw;
        
        conudctHeat = kc*(T(2,i)-T(1,i))/(H(i)*dzeta);
        m_basal(i) = (Qgeo + conudctHeat + frictionHeat(i))/(rhow*Lw)*SPY; % basal melt [m a-1]
                
        % CTS jump condition
        if jcts>1 & w_vs(jcts,i)<0
            E(jcts,i) = Epmp_i(jcts);
        end
    end
        
    % diffusive water flux through temperate ice
    qw_TEMP_diffu(:,i) = -Kt*(omega(2:end,i)-omega(1:end-1,i))/(H(i)*dzeta);
    
    % vertically integrated drainage [m s-1]
    % drain instantaneously to the bed using the Ralf Greve's function
    qw_greveDrain(i) = sum(drainageRate)*H(i)*dzeta; % [m s-1]
    
    % set the upper bound of the water content (3%)
    %     index1 = omega(:,i) > 3/100;
    %     index2 = omega(:,i) <= 3/100;
    %     omega(:,i) = index1*3/100 + index2.*omega(:,i);
end

end