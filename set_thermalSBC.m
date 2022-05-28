function [Esbc] = set_thermalSBC(para)
% set the surface boundary condition for the thermal model
% 2019-2-15

global hS M

Cp = para.Cp;
Tref = para.Tref;

%% Storglaciaeren example
%
% T0 = 273.15;
% Tma =  -6.0; % degC, mean annual air temperature at Tarfala
% zcts = 1300; % m a.s.l.; altitude where CTS is at the surface, projected to topg
% slope = 100; % m; range around which surface temp transition happens
% Tsbc = T0 + Tma * (zcts + slope - hS) / (2.0 * slope);
% Tsbc(hS<zcts-slope) = T0 + Tma;
% Tsbc(hS>zcts+slope) = T0;

%% Geladandong
% 
% T0 = 273.15-4.5;
% z0 = 5797;
% 
% z_ela = median(hS);
% Tsbc = zeros(1,M);
% for i=1:M
%     if hS(i)<=z_ela
%         Tsbc(i) = T0 - 8.75e-3*(hS(i)-z0);
%     else
%         Tsbc(i) = 273.15 - 6.5;
%     end
% end

%% temperate ice experiment

T0 = 273.15 - 20.0;
z0 = hS(end);

z_ela = median(hS);
Tsbc = zeros(1,M);
for i=1:M
    if hS(i)<=z_ela
        Tsbc(i) = T0 - 6.5e-3*(hS(i)-z0);
    else
        Tsbc(i) = T0 - 6.5e-3*(hS(i)-z0);
    end
end

% Tsbc = (273.15 - 1)*ones(1,M);
%% Enthalpy
Esbc = Cp*(Tsbc - Tref);
end
