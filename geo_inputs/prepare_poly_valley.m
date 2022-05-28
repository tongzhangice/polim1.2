clear,clc

load geo_arolla.mat
hS = geo.hS;
hB = geo.hB;
xi = geo.xi;
Wsuf = geo.Wsurf;

delta_hS = max(hS) - min(hS);
delta_hB = max(hB) - min(hB);
delta_xi = max(xi) - min(xi);

norm_hS = (hS-min(hS))/delta_hS;
norm_hB = (hB-min(hB))/delta_hB;
norm_xi = (xi-min(xi))/delta_xi;

hS1 = norm_hS*1200 + 4000;
hB1 = norm_hB*1200 + 4000;
xi1 = norm_xi*10e3;
H1 = hS1 - hB1;

hold on
plot(xi, hS, 'k-', 'linewidth', 1)
plot(xi, hB, 'k-', 'linewidth', 1)

plot(xi1, hS1, 'r-', 'linewidth', 1)
plot(xi1, hB1, 'r-', 'linewidth', 1)

geo.xi = xi1;
geo.hS = hS1;
geo.hB = hB1;
geo.H = H1;

% save geo_poly_valley geo
