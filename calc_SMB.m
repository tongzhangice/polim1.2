% Date: 2020/11/29

function SMB_ice = calc_SMB(zref, Tma, hSn)
M = length(hSn);
rhow = 1000.0;
rhoi = 910.0;

d = 1:365; % days
d = d';
%% Surface air temperature
%
% 
fdd = 7.0e-3; % degree-day factor [m K-1 d-1]
lapseRate = -0.0075; % lapse rate [K m-1]

deltaT = 0; % -4 ~ 4 K

% daily air temperature at the reference elevation
Tref = -10*cos(2*pi*d/365) + Tma + deltaT;

% Surface air temperature along the flowline
Tair = Tref*ones(1,M) + ones(365,1)*(hSn-zref)*lapseRate;

dailyMelt = fdd*Tair; % [m s-1]
dailyMelt(dailyMelt<0) = 0;

%% Precipitation
%
% 
dmax = 168; % June 16
Ia = 12; % parameter about seasonality concentration []
Iw = 10; % parameter about weekly amplitude []
grad_P = 2.2e-4; % precipitation gradient [m / m]

% precipitation ratio
rario_P = Ia*cos((d-dmax)*2*pi/365) + Iw*cos((d-dmax)*2*pi*52/365);
rario_P = max(rario_P, 0);
sum_Pr = sum(rario_P);

% daily precipitation
annualPrep = (hSn-zref)*grad_P + 0.3; % annual mean precipitation along the flowline
dailyPrep = rario_P/sum_Pr*annualPrep;

%% Surface mass balance
%
%
SMB = sum(dailyPrep) - sum(dailyMelt);
SMB_ice = SMB*rhow/rhoi;

end
