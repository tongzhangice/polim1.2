clear,clc
%% ISMIP-HOM Exp.D, ice stream
% iceL = 10e3;
% geo.xi = linspace(0,iceL,51);
% slp = 0.1;
% geo.hS = -geo.xi.*tan(slp*pi/180);
% geo.hB = geo.hS - 1e3;
% geo.H = geo.hS - geo.hB;
% 
% save geo_ismip_expD_L geo


%% Polythermal ice slab
% geo.xi = 100;
% geo.dx = 100;
% geo.hS = 200;
% geo.hB = 0;
% geo.H = 200;
% save geo_polythermal_slab geo

%% valley glacier (Hewitt&Schoof, 2017)
% ice_L = 10e3;
% dx = 1e3;
% xi = linspace(0,ice_L, 51);
% hB = 800*(1-xi/ice_L) - 800*(1-xi/ice_L).*(xi/ice_L) + 4000;
% hS = 800*(1-xi/ice_L) + 300*(1-xi/ice_L).*(xi/ice_L) + 4000;
% H = hS - hB;
% geo.xi = xi;
% geo.hS = hS;
% geo.hB = hB;
% geo.H = H;
% save geo_hewitt_valley geo
% 
% figure(1);
% plot(xi, hS, 'b-','linewidth',2)
% hold on
% plot(xi, hB, 'k-', 'linewidth',2)
% figure(2); plot(xi, H, 'k-')

% extend valley glacier (Hewitt&Schoof, 2017)
% ice_L = 6e3;
% dx = 1e2;
% xi = 0:dx:ice_L;
% hB = 500*(1-xi/ice_L) - 500*(1-xi/ice_L).*(xi/ice_L) + 4000;
% hS = 500*(1-xi/ice_L) + 300*(1-xi/ice_L).*(xi/ice_L) + 4000;
% 
% Hmin = 0.5;
% hB1 = [hB, hB(end)*ones(1,8)];
% hS1 = [hS, hB(end)*ones(1,8)+Hmin];
% xi1 = [xi, (xi(end)+dx):dx:(xi(end)+8*dx)];
% H1 = hS1 - hB1;
% 
% geo.xi = xi1;
% geo.hS = hS1;
% geo.hB = hB1;
% geo.H = H1;
% 
% save geo_hewitt_valley_extend geo
% figure(1); plot(xi1, hS1, xi1, hB1)

%% Ice cap (Hewitt&Schoof, 2017)
% ice_L = 100e3;
% dx = 1e3;
% xi = 0:dx:ice_L;
% hB = zeros(1,length(xi));
% hS = 1500.*(1 - (xi/ice_L).^2);
% H = hS - hB;
% geo.xi = xi;
% geo.hS = hS;
% geo.hB = hB;
% geo.H = H;
% save geo_hewitt_icecap geo
%% Slab for Robin inversion
% dx = 100;
% xi = 0:dx:5000;
% hB = 100 - 0.1*xi;
% hS = hB + 200;
% H = hS - hB;
% 
% geo.xi = xi;
% geo.hS = hS;
% geo.hB = hB;
% geo.H = H;
% save geo_slab geo

%% Beaud (2014)
% ice_L = 10e3;
% dx = 200;
% xi = 0:dx:ice_L;
% 
% hS = -400/ice_L^2*xi.^2 + 400;
% hB = zeros(1,length(xi));
% H = hS - hB;
% 
% geo.xi = xi;
% geo.hS = hS;
% geo.hB = hB;
% geo.H = H;
% 
% % plot(geo.xi, geo.hS, geo.xi, geo.hB)
% save geo_beaud geo

%% Hoffman (2014)
% rhoi = 910;
% g = 9.81;
% tau_c = 1e5;
% ice_L = 7e3;
% b0 = 700;
% 
% xi1 = linspace(ice_L, 0, 701);
% hB = b0*(1-xi1/ice_L);
% 
% xspan = [ice_L, 0];
% opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
% hS0 = 2;
% [x_solve,hS] = ode45(@(x,hS) ode_NyePlastic(x,hS,xi1,hB), xspan, hS0, opts);
% 
% xi = linspace(0, ice_L, 51);
% hS = interp1(x_solve, hS, xi);
% hB = interp1(xi1, hB, xi);
% 
% H = hS - hB;
% 
% geo.xi = xi;
% geo.hS = hS;
% geo.hB = hB;
% geo.H = H;
% save geo_nye geo
% 
% figure
% plot(xi, hS, 'b-')
% hold on
% plot(xi, hB, 'k-')
% 
% xlabel('x (km)')
% ylabel('Elevation (m)')
% set(gcf, 'position', [430.0783  379.1565  651.1304  243.7565])

%% SHMIP Icesheet
ice_L = 100e3;
geo.xi = linspace(0,ice_L,51);
dx = geo.xi(2)-geo.xi(1);

hS = 6*((geo.xi+5000).^(1/2)-sqrt(5000))+1;
hB = zeros(1,length(geo.xi));
geo.hS = fliplr(hS);
geo.hB = fliplr(hB);
geo.H = geo.hS - geo.hB;
geo.Wsurf = 10000*ones(1,length(geo.xi));
save geo_shmip_icesheet geo

%% SHMIP Glacier
% ice_L = 6e3;
% dx = 100;
% geo.xi = 0:dx:ice_L;
% M = length(geo.xi);
% 
% gammab = 0.05;
% hS = 100*(geo.xi+200).^(1/4) + geo.xi/60 - (2*10^10)^(1/4) + 1;
% hB = (hS(end) - 6000*gammab) / (6000^2) * geo.xi.^2 + gammab * geo.xi;
% geo.hS = fliplr(hS);
% geo.hB = fliplr(hB);
% geo.H = geo.hS - geo.hB + 0.1;
% geo.Wsurf = 10000*ones(1,M);
% save geo_shmip_glacier geo

%% ISMIP-HOM Haut d'Arolla
% geoGlacier = dlmread('arolla100.dat');
% geoGlacier = geoGlacier';
% 
% % origin
% geo.xi = geoGlacier(1,:);
% geo.hB = geoGlacier(2,:);
% geo.hS = geoGlacier(3,:);
% geo.H = geo.hS - geo.hB + 0.1;

%% ISMIP-HOM Haut d'Arolla resample
% geoGlacier = dlmread('arolla100.dat');
% geoGlacier = geoGlacier';        
% % origin
% xi0 = geoGlacier(1,:);
% hB0 = geoGlacier(2,:);
% hS0 = geoGlacier(3,:);
% 
% dx = 10;
% xi = 0:dx:5000;
% hB1 = interp1(xi0, hB0, xi, 'spline');
% hS1 = interp1(xi0, hS0, xi, 'spline');
% 
% hB2 = smooth(hB1, 0.1, 'loess');
% hS2 = smooth(hS1, 0.1, 'loess');
% 
% hS = hS2'; hB = hB2';
% 
% logic1 = (hS-hB) > 0;
% hB = logic1.*hB + (1-logic1).*(hS - 0.1);
% H = hS - hB;
% 
% geo_arolla_resample.xi = xi;
% geo_arolla_resample.hS = hS;
% geo_arolla_resample.hB = hB;
% geo_arolla_resample.H = H;
% 
% save geo_arolla_resample geo_arolla_resample