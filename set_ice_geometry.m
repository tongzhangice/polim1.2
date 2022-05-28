function set_ice_geometry(geofile, para)

global xi dx hS hB H W W_s M Ms N dzeta zeta

type_valley = para.type_valley;

%% LOAD
%
%
load(geofile)
xi = geo.xi;
hB = geo.hB;
hS = geo.hS;
H = geo.H + para.Hmin;
% H = geo.H;
% H(end-1:end) = H(end-1:end) + para.Hmin;

%
% MODEL GRID
if size(xi)==1
    dx = 100;
else
    dx = xi(2) - xi(1);
end
M = length(xi);
Ms = M + 1;

N = para.layers;
dzeta = 1/(N-1);
zeta = 0:dzeta:1;
zeta = zeta';

%% WIDTH
%
Wsurf_const = 2e3;

% realistic plan-view half flowband width
if isfield(geo, 'Wsurf')
    if size(geo.xi)==size(geo.Wsurf)
        Wsurf = geo.Wsurf;
    else
        Wsurf = Wsurf_const/2*ones(1,M);
    end
else
    Wsurf = Wsurf_const/2*ones(1,M);
end

% rectangular plan-view
% Hewitt valley, temperate ice exp, Wsurf=2e3
% ismip expD inverse, Wsurf=50e3
% Wsurf = 20e3/2*ones(1,M);

% trapezoid plan-view
% Whead = 500;
% Wterm = 50;
% Wsurf = (Whead-Wterm)*(xi(end)-xi)/xi(end) + Wterm;

%
% Half-width distribution
W = zeros(N,M);
if strcmpi(type_valley,'sves')
    para_sves = 2*ones(1,M);
    for i=1:M
        W(:,i) = Wsurf(i)*(zeta.^(1/para_sves(i)));
    end
    W(1,:) = W(2,:)/2;
elseif strcmpi(type_valley, 'rect')
    W = ones(N,1)*Wsurf;
elseif strcmpi(type_valley,'trapz')
    alpha_trapezoid = 45*pi/180;
    W = ones(N,1)*Wsurf - (ones(N,1)*H - zeta*H)/tan(alpha_trapezoid);
    W(W<1e-3) = 1;
else
    disp('p.type_valley should be svesson, rect, trapz!')
end

W_s = main2staggerX(W);

end