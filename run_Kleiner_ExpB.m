clc
clearvars
clearvars -global

tic
global N M H dzetadx zeta iTimeStep...
    At_E At_T At_omega At_CTS At_Kappa_vs...
    At_isTEMP At_thk_TEMP At_thk_w    

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params_Kleiner_ExpB();
SPY = p.SPY;

p.layers = 21;
%% GLACIER GEOMETRY
%
%
set_ice_geometry('./geo_inputs/geo_Kleiner_ExpB.mat', p);
dzetadx = zeros(N,M);
%% TIME SETTING
%
%
dt_u = 1; % time step for velocity solver [a]
endTime = 2000; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
alpha = 4*pi/180; % inclination angle
vz = -0.2; % [m a-1] vertical velocity

u = zeros(N,M);
u_s = zeros(N,M+1);
w = vz*ones(N,M);
w_vs = vz*ones(N,M);
strainHeat = 2*p.AGlen*(p.rhoi*p.g*sin(alpha)).^4.*H.^4.*(1-zeta).^4*SPY;
frictionHeat = zeros(1,M);

%  Related to enthalpy solver
At_E = zeros(N,M,numTimeStep);
At_T = zeros(N,M,numTimeStep);
At_omega = zeros(N,M,numTimeStep);
At_Kappa_vs = zeros(N-1,M,numTimeStep);
At_CTS = zeros(numTimeStep,M);
At_isTEMP = zeros(numTimeStep,M);
At_thk_TEMP = zeros(numTimeStep,M);
At_thk_w = zeros(numTimeStep,M); % water thickness
At_m_basal = zeros(numTimeStep,M); % basal melt rate

% SBC and IC for the enthalpy solver
Esbc = p.Cp*(-3 + 273.15 - p.Tref);
Eini = p.Cp*(-1.5 + 273.15 - p.Tref)*ones(N,M);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    trueTime = arrayTime(iTimeStep);
    
    % calculate the enthalpy
    [E, T, omega, Kappa_vs, CTS, thk_TEMP, m_basal, qw_greveDrain, qw_TEMP_diffu] = ...
        solver_enthalpy_SEGM(u, u_s, w, w_vs, strainHeat, frictionHeat, dt_u, Esbc, Eini, p);
    
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
end
toc

% figure
plot_enthalpy_ExpB