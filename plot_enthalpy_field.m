global xi M N zeta dzeta H hS hB
SPY = 31536000;

% istep = 800;

if enth_nan
    istep = iTimeStep-1;
else
    istep = iTimeStep;
end
    
xx = ones(N,1)*xi;
yy = ones(N,1)*hB + zeta*H;


if length(size(At_CTS))==2
    CTS_istep = At_CTS(istep,:);
else
    CTS_istep = At_CTS;
end

if length(size(At_E))==3  
    E_istep = At_E(:,:,istep);
    T_istep = At_T(:,:,istep);
    omega_istep = At_omega(:,:,istep);
    u_istep = At_u(:,:,istep);
else
    E_istep = At_E;
    T_istep = At_T;
    omega_istep = At_omega;
    u_istep = At_u;
end

index1 = find(CTS_istep>0 & CTS_istep<N);

if ~isempty(index1)
    if index1(1)==1
        if index1(end)==M
            index2 = index1; % 1:M
        else
            index2 = [index1, index1(end)+1];
        end
    else
        if index1(end)==M
            index2 = [index1(1)-1, index1];
        else
            index2 = [index1(1)-1, index1, index1(end)+1];
        end
    end
end
%% figure
%
figure

%% Enthalpy field
%
ax(1) = subplot(2,2,1);
hold on
contourf(xx/1000, yy, E_istep, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
    plot(xi(index2)/1000, hB(index2)+CTS_istep(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(1), flipud(hot))

ylabel(c, 'Enthalpy (J kg^{-1})')

xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3])
box on

%% Porosity field
%
ax(2) = subplot(2,2,2);
hold on
contourf(xx/1000, yy, omega_istep*100, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
plot(xi(index2)/1000, hB(index2)+CTS_istep(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(2), brewermap([],'*RdYlBu')) % *RdYlBu *RdBu
% colormap(ax(3), 'jet')

ylabel(c, 'Moisture (%)')

xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3])

box on

%% Velocity
%
ax(3) = subplot(2,2,3);
hold on
contourf(xx/1000, yy, u_istep, 50, 'LineStyle', 'none')
plot(xi/1000, hB, 'k', 'linewidth', 2)
plot(xi/1000, hS, 'k', 'linewidth', 2)

if ~isempty(index1)
plot(xi(index2)/1000, hB(index2)+CTS_istep(index2).*H(index2)*dzeta, 'w-', 'linewidth',1.5)
end
hold off

c = colorbar;
colormap(ax(3), 'jet')

ylabel(c, 'u (m a^{-1})')

xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
xlim([0, xi(end)/1e3])

box on

%% Darcy water flux and basal melt
%
ax(4) = subplot(2,2,4);
plot(xi/1e3, At_qw_TEMP_darcy(1,:,istep)*SPY*1e3, 'k-')
hold on
plot(xi/1e3, -m_basal*1e3, 'r-')
xlabel('Horizontal distance (km)', 'FontSize', 10)
ylabel('Water flux (mm a^{-1})', 'FontSize', 10)
box on
grid on

legend('Darcy flux', 'Basal melt', 'Location', 'SouthWest')
set(gcf, 'Position', [558,275,831,467])


% ax(5) = subplot(3,2,5);
% nyear = size(At_CTS,1);
% arr_color = jet(nyear);
% hold on
% if nyear==1
%     plot(xi/1000, hB+At_CTS(1,:).*H*dzeta, 'Color', 'r')
% else
%     for i=1:nyear
%         plot(xi/1000, hB+At_CTS(i,:).*H*dzeta, 'Color', arr_color(i,:))
%     end
%     c = colorbar;
%     colormap(arr_color)
%     caxis([1,nyear])
%     ylabel(c, 'Year')
% end
% 
% plot(xi/1000, hB, 'k', 'linewidth', 2)
% % plot(xi/1000, hS, 'k', 'linewidth', 2)
% hold off
% box on
% xlabel('Horizontal distance (km)', 'FontSize', 10)
% ylabel('Elevation (m a.s.l.)', 'FontSize', 10)
% xlim([0, xi(end)/1e3])
% % axis tight

% set(gcf, 'Position', [362,151,796,602])