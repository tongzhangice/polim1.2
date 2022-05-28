clc
clearvars
clearvars -global

tic
global N M hB H iTimeStep zeta dhSdx...
    At_E At_T At_omega At_CTS At_Kappa_vs...
    At_isTEMP At_thk_TEMP At_thk_w At_qw_TEMP_darcy

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params_hewitt_icecap();
SPY = p.SPY;
p.Hmin = 0;

p.type_enth = 'MEGM';
p.is_auto_enth_BBC = 0; % At bed: T = Tm
p.type_enth_BBC = 2;
p.has_Greve_drainage = 1;

p.layers = 101;

%% GLACIER GEOMETRY
%
%
set_ice_geometry('./inputs/geo_hewitt_icecap', p);

%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 4000; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
%-----------------solve the velocity field-----------------
set_staggered_grid();

u = zeros(N,M);
strainHeat = zeros(N,M);
for i = 1:M
    for j = 1:N
        u(j,i) = -2*p.AGlen*(p.rhoi*p.g*dhSdx(i)).^3/4*(H(i).^4 - ((1-zeta(j))*H(i)-hB(i)).^4);
        strainHeat(j,i) = 2*p.AGlen*(p.rhoi*p.g*dhSdx(i)).^4*((1-zeta(j)).*H(i)-hB(i)).^4;
    end
end
u_s = main2staggerX(u);
[w_vs, w] = get_ice_w(u_s, u);

% convert [m s-1] to [m a-1]
u = u*SPY;
u_s = u_s*SPY;
w = w*SPY;
w_vs = w_vs*SPY;
strainHeat = strainHeat*SPY;
%-----------------solve the velocity field-----------------

%  Related to enthalpy solver
At_E = zeros(N,M,numTimeStep);
At_T = zeros(N,M,numTimeStep);
At_omega = zeros(N,M,numTimeStep);
At_Kappa_vs = zeros(N-1,M,numTimeStep);
At_CTS = zeros(numTimeStep,M);
At_thk_TEMP = zeros(numTimeStep,M);
At_m_basal = zeros(numTimeStep,M); % basal melt rate
At_thk_w = zeros(numTimeStep,M); % water thickness
At_Tb = zeros(numTimeStep,M); % basal temperature
At_isTEMP = zeros(numTimeStep,M);
At_qw_TEMP_darcy = zeros(N-1,M,numTimeStep);

enth_covg = 0;
enth_nan = 0;
pdtj_storage = [];

% Related to glacier evolution
At_hS = zeros(numTimeStep,M);
At_hB = zeros(numTimeStep,M);
At_H = zeros(numTimeStep,M);
At_SMB = zeros(numTimeStep,M);

% SBC and IC for the enthalpy solver
Ts = -1;
Esbc = p.Cp*(Ts + 273.15 - p.Tref)*ones(1,M);
Eini = p.Cp*(262.15 - p.Tref)*ones(N,M);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    % calculate the enthalpy
    if strcmpi(p.type_enth,'SEGM')
        [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_greveDrain, qw_TEMP_diffu] = ...
            solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth,'MEGM')
        [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_TEMP, qw_TEMP_darcy] =...
            solver_enthalpy_MEGM(u, u_s, w, w_vs, strainHeat, dt_u, Esbc, Eini, p);
    elseif strcmpi(p.type_enth,'isoT')
        % pass
    end

    % converge
%     if iTimeStep>1
%         Elst = At_E(:,:,iTimeStep-1);
%         pdtj = sumsqr(E-Elst)/sumsqr(Elst);
%         pdtj_storage = [pdtj_storage; pdtj];
%         if pdtj<2.1e-8
%             enth_covg = 1;
%             break
%         end
%     end
    
    %%
    %
    % STORAGE RESULTS
    At_E(:,:,iTimeStep) = E;
    At_T(:,:,iTimeStep) = T;
    At_omega(:,:,iTimeStep) = omega;
    At_Kappa_vs(:,:,iTimeStep) = Kappa_vs;
    At_CTS(iTimeStep,:) = CTS;
    At_thk_TEMP(iTimeStep,:) = thk_TEMP;
    At_m_basal(iTimeStep,:) = m_basal;
    
    % calculate the basal water layer thickness
    if iTimeStep == 1
        At_thk_w(iTimeStep,:) = zeros(1,M) + dt_u*m_basal;
    else
        At_thk_w(iTimeStep,:) = At_thk_w(iTimeStep-1,:) + dt_u*m_basal;
    end
    
    % logic value for basal thermal state
    logic1 = (At_thk_w(iTimeStep,:)>0 & At_isTEMP(iTimeStep,:)==0);
    At_isTEMP(iTimeStep,:) = At_isTEMP(iTimeStep,:) | logic1;
    
    if strcmpi(p.type_enth,'MEGM')
        At_qw_TEMP_darcy(:,:,iTimeStep) = qw_TEMP_darcy;
    end
end
toc

plot_hewitt_icecap