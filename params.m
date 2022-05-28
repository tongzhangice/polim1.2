function p = params()

% PHYSICAL CONSTANTS
p.g = 9.81; % accel of gravity [m s-2]
p.R = 8.314; % universal gas constant [J mol K-1];
p.SPD = 86400; % day is this many seconds
p.SPY = 31536000; % 31556926; % year is this many seconds (i.e. 365.2422 days)
p.rhoi = 910.0; % density of ice [kg m-3]
p.rhow = 1000.0; % water density [kg m-3]
p.n = 3; % exponent in Glen's flow law []
p.de0 = 1e-30; % small number in case of singularity [yr-2]
p.Lw = 3.34e5; % latent heat of fusion [J kg-1]
p.betaCC = 7.9e-8; % Clausius-Clapeyron constant [K Pa-1]
p.eta_w = 1e-3; % water viscosity (Hewitt&Schoof, 2017, 1.8e-3) [Pa s]
p.iter_max = 50; % maximum iterations for velocity solver
p.Hmin = 0.1; % minimal ice thickness in case of singularity
p.f_nye = 1; % Nye shape factor
p.layers = 31; % vertical layers

% parameters related to thermodynamics
p.Tref = 223.15; % Reference temperature [K]
p.kc = 2.1; % cold ice conductivity [W m-1 K-1]
p.Cp = 2009; % ice specific heat capacity [J kg-1 K-1]
p.Kc = p.kc/(p.rhoi*p.Cp); % thermal diffusivity of cold ice [m2 s-1]
p.Kt = 1.1e-11; % regularization of thermal diffusivity of temperate ice [m2 s-1]
p.k0 = 1e-12; % permeability factor (unconstrained, Hewitt&Schoof, 2017, Tab. 1) [m2]; range: 1e-12~5e-8
p.alpha_hewitt = 2; % exponent in compaction pressure model (Hewitt&Schoof, 2017); unconstrained, range: 2~3
p.Qgeo = 0.06; % geothermal heat flux [W m-2]
p.moi_UB = 0.01; % upper bound of moisture in temperate ice

% parameters related to Coulomb sliding law
p.lambda_max = 4; % wavelength of the dominant bedrock bumps [m]
p.m_max = 0.5; % maximum slope of the dominant bedrock bumps []; ref: 0.5
p.kflot = 0.1; % a fraction of flotation; N=Pi-Pw, Pw=kflot*Pi, N=(1-kflot)*Pi

% parameters related to subglacial hydrology
p.lr = 2; % wavelength of bedrock bumps [m]; ref: 2
p.hr = 0.1; % height of bedrock bumps [m]; ref: 0.1
p.ktransm = 4e-4; % transmissivity coefficient [], larger value, larger N; ref: 4e-4
p.ev = 1e-5; % englacial void fraction (de Fleurian, 2018, eq. 10)
p.Hw_crit = 0.15; % critcal water sheet thickness (Pimentel & Flowers, 2010, Tab. 1)

% parameters related to surface mass balance
p.SMB_grad = 4e-3;
p.SMB_max = 1.5;

% OPTIONS
p.is_iceflow = 1;

p.is_flowband = 1;
% true  (1): flowband mode
% false (0): flowline mode

p.is_thk_evolv = 0;
% true (1): ice thickness evolution in transient simulation
% false (0): ice thickness is fixed

p.is_surf_relax = 0; % surface relaxation

p.is_subHyro = 1;
% true (1): couple the subglacial hydrology to ice dynamics
% false (0): without considering the subglacial hydrology

p.is_auto_enth_BBC = 1;
% true (1): basal thermal BC is automatically set depending on Aschwanden
% decision chart (Aschwanden et al., 2012, JG, p451, fig. 5)
% false (0): basal thermal BC is specified

p.is_duval = 1;
% true (1): considering the effect of water on temperate ice creep based on Duval-Lliboutry relation
% false (0): without considering the effect of water on temperate ice creep

p.type_SBC = 'neum';
% neumann: stress-free boundary condition
% dirichlet: specified surface velocity

p.type_BBC = 'LFlaw_simple'; % ;
% zero: no-slip bed
% CFlaw_polyT: Coulomb friction law, sliding depends on basal thermal state
% CFlaw_isoT: Coulomb friction law, sliding everywhere
% LFlaw: linear friction law
% LFlaw_E2: linear friction law for ISMIP-HOM E2
% LFlaw_simple: simplified linear friction law

p.type_valley = 'rect';
% trapz: trapezoid
% sves: Svessen, y=ax^b
% rect: rectangular

p.type_LBC = 'zero';
% zero: zero velocity at the terminus
% calv: calving front

p.type_enth_solver = 'isoT';
% SEGM: standard enthalpy gradient model (SEGM; Aschwanden et al., 2012)
% MEGM: modified enthalpy gradient model (MEGM; Hewitt I. and Schoof C., 2017)
% isoT: isothermal

p.type_enth_BBC = 1;
% see Kleiner et al., 2015, TC, section 2.2
% 1: Cold base (dry)
% 2: Temperate base
% 3: Temperate ice at base
% 4: Cold base (wet)

p.has_Greve_drainage = 1;
% This option is only used for SEGM model.
% see Aschwanden et al., 2012, JG, p450, fig. 4
% 1: use Greve drainage function
% 0: do not use Greve drainage function

p.type_Arrhenius = 'cuffey_duvalON';
% greve: Greve (2009)
% cuffey: Cuffey & Paterson (2010)

p.is_SBC_width = 0;
% true  (1): incorporate flowband parameterization in stress-free condition
% false (0): do not considier flowband parameterization in stress-free condition
end