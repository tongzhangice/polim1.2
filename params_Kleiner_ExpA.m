% 2019-1-16 
function p = params_Kleiner_ExpA()

p.SPY    = 31556926;    % seconds in one year [s a-1]
p.rhoi   = 910;         % ice density [kg m-3]
p.rhow   = 1000;        % water density [kg m-3]
p.g      = 9.81;        % gravitational acceleration [m s-2]
p.kc     = 2.1;         % cold ice conductivity [W m-1 K-1]
p.Cp     = 2009;        % ice specific heat capacity [J kg-1 K-1]
p.Lw     = 3.34e5;      % latent heat of fusion [J kg-1]
p.Tref   = 223.15;      % Reference temperature [K] 273, 223.15;
p.Kc     = p.kc/(p.rhoi*p.Cp); % Thermal diffusivity of cold ice [m2 s-1]
p.Kt     = p.kc/(p.rhoi*p.Cp)*1e-1;      % Thermal diffusivity of temperate ice [m2 s-1]
p.betaCC = 7.9e-8;           % Clausius-Clapeyron constant [K Pa-1]
p.Qgeo   = 0.042;
p.n      = 3;

p.type_enth_solver = 'SEGM';
p.type_valley = 'rect';
p.Hmin = 0.0;
p.has_Greve_drainage = 0;
p.is_auto_enth_BBC = 1;
p.type_enth_BBC = 1;
p.layers = 201;
end