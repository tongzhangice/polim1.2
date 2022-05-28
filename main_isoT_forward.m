clc
clearvars
clearvars -global

tic
global N M Ms hS hB H iter_u u_s_lst iTimeStep At_u xi

%% PHYSICAL CONSTANTS AND PARAMETERS
%
%
p = params();
SPD = p.SPD;
SPY = p.SPY;
iter_max = p.iter_max;

p.is_flowband = 1;
p.type_BBC = 'LFlaw'; %'LFlaw_simple'; %'zero'; %'CFlaw_isoT';
p.layers = 31;
p.Hmin = 0.1;

p.type_valley = 'rect';

%% GLACIER GEOMETRY
%
%
set_ice_geometry('./inputs/geo_arolla', p); % geo_shmip_icesheet, geo_beaud



%% TIME SETTING
%
%
dt_u = 1*SPY; % time step for velocity solver [a]
endTime = 1*SPY; % end time of model run [a]
[arrayTime, numTimeStep] = set_time_step(dt_u, endTime);

%% INITIALIZATION
%
%
AGlen_s = zeros(N,Ms) + 1e-16; % [Pa-3 a-1]
visc_s = zeros(N,Ms) + 1e13/SPY; % [Pa a]
visc = zeros(N,M) + 1e13/SPY; % [Pa a]

%  Related to velocity solver
At_u = zeros(N,M,numTimeStep);
At_w = zeros(N,M,numTimeStep);

%% MAIN
%
%
for iTimeStep = 1:numTimeStep
    fprintf('iTimeStep: %d \n', iTimeStep)
    
    
    set_staggered_grid();
    
    iter_u = 0;
    while 1
        iter_u = iter_u + 1;
        fprintf('iter_u: %d\n', iter_u)
        
        [u_s] = solver_u(visc_s, visc, AGlen_s, p);
        %----------------------Begin Picard iteration----------------------
        if iter_u>2
            u_s_now = u_s;
            Cs = u_s_now - u_s_lst;
            Sita = acos(Cs'*C/(sumsqr(Cs)*sumsqr(C)));
            if isequal(Sita<=pi/8, ones(Ms,Ms))
                mu1 = 2.5;
            elseif isequal(Sita>pi/8,ones(Ms,Ms)) && isequal(Sita<19*pi/20, ones(Ms,Ms))
                mu1 = 1;
            elseif isequal(Sita>=19*pi/20, ones(Ms,Ms))
                mu1 = 0.5;
            end
            u_s = u_s_lst + mu1*Cs;
            
            if sumsqr(u_s_now - u_s_lst)/sumsqr(u_s_now)<1e-4
                break
            end
        end
        
        if iter_u>1
            C = u_s - u_s_lst;
        end
        
        if iter_u>=iter_max
            break
        end
        
        u_s_lst = u_s;
        %-----------------------End Picard iteration-----------------------
        u = staggerX2main(u_s);
        
        % calculate the vertical velocity
        [w_vs, w] = get_ice_w(u_s, u);
        
        % calculate the strain heat
        [visc_s, visc, strainHeat] = get_ice_viscosity(u_s, u, AGlen_s, p);
    end
    % calculate basal shear stress and basal melting rate
    [tauxz, tauxz_s] = calc_tauxz(u_s, visc_s); % [Pa]
    %         m_basal = (p.Qgeo + tauxz.*u(1,:)/SPY)/(p.rhoi*p.Lw); % [m s-1]
    
    fprintf('Max sliding velocity: %3.2f \n', max(u(1,:)))
    fprintf('Max surface velocity: %3.2f \n', max(u(end,:)))
    At_u(:,:,iTimeStep) = u;
    At_w(:,:,iTimeStep) = w;
end

toc

plot_uField_uSurf