global xi N H zeta hS hB

% if enth_nan
%     istep = iTimeStep-1;
% else
%     istep = iTimeStep;
% end

istep = iTimeStep;

figure
subplot(2,1,1)
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;

if length(size(At_u))==3  
    u_istep = At_u(:,:,istep);
else
    u_istep = At_u;
end
    

contourf(xx/1000, yy, u_istep, 50, 'LineStyle', 'none')
hold on
plot(xi/1000, hB, 'k')
plot(xi/1000, hS, 'k')

vcb = colorbar;
set(vcb, 'location', 'North', 'position', [0.63 0.84 0.25 0.02], 'fontsize', 9)
colormap('Jet')

ylabel(vcb, 'Velocity (m a^{-1})', 'fontsize', 9)

xlabel('Horizontal distance (km)')
ylabel('Elevation (m a.s.l.)')

% load('geo_arolla.mat')
% plot(geo.xi/1000, geo.hS, 'k-', 'linewidth', 1.5)
% plot(geo.xi/1000, 3940*ones(1,M), 'k--', 'linewidth', 2)

subplot(2,1,2)
HL = zeros(1,2);
HL(1) = plot(xi/1000, u_istep(end, :), 'k-', 'linewidth', 3);
hold on
HL(2) = plot(xi/1000, u_istep(1,:), 'r-', 'linewidth', 0.5,...
    'Marker', '.', 'MarkerSize', 10);

xlabel('Horizontal distance (km)')
ylabel('Velocity (m a^{-1})')

hlgd = legend(HL, 'surface velocity', 'sliding velocity');
set(hlgd, 'fontsize', 8, 'location', 'NorthWest');

% ylim([-3, max(max(u_istep))+20])
% title('Surface velocity')

set(gcf, 'Position', [661,203,542,615])

% axis tight