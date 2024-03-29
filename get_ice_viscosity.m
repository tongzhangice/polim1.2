function [visc_s, visc, strain_heat] = get_ice_viscosity(u_s, u, AGlen_s, para)
% Inputs:
    % u_s [m a-1]
    % u [m a-1]
    % AGlen_s [Pa-3 a-1]
% Outputs:
    % visc_s [Pa a]
    % visc [Pa a]
    % strain_heat [Pa a-1]
% global isFlowband de0 n
global Ms N dzeta dx dzetadx_s H_s W W_s

is_flowband = para.is_flowband;
de0 = para.de0;
n = para.n;

de2_s=zeros(N,Ms); % [a-2]
for i = 1:Ms
    for j = 2:N-1        
        if i==1
            de2_s(j,i) = ((u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-8*u_s(j,i)+9*u_s(j,i+1)-u_s(j,i+2))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2 + ...
                (u_s(j,i)/W_s(j,i))^2*((-8*W_s(j,i)+9*W_s(j,i+1)-W_s(j,i+2))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta))^2 + ...
                u_s(j,i)/W_s(j,i)*((-8*u_s(j,i)+9*u_s(j,i+1)-u_s(j,i+2))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta)) *...
                ((-8*W_s(j,i)+9*W_s(j,i+1)-W_s(j,i+2))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(j,i)/W_s(j,i))^2;
            
            de2_s(N,i) = ((3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-8*u_s(N,i)+9*u_s(N,i+1)-u_s(N,i+2))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2 + ...
                (u_s(N,i)/W_s(N,i))^2*((-8*W_s(N,i)+9*W_s(N,i+1)-W_s(N,i+2))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta))^2 + ...
                u_s(N,i)/W_s(N,i)*((-8*u_s(N,i)+9*u_s(N,i+1)-u_s(N,i+2))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta)) * ...
                ((-8*W_s(N,i)+9*W_s(N,i+1)-W_s(N,i+2))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(N,i)/W_s(N,i))^2;
            
            de2_s(1,i) = ((-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-8*u_s(1,i)+9*u_s(1,i+1)-u_s(1,i+2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2 + ...
                (u_s(1,i)/W_s(1,i))^2*((-8*W_s(1,i)+9*W_s(1,i+1)-W_s(1,i+2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta))^2 + ...
                u_s(1,i)/W_s(1,i)*((-8*u_s(1,i)+9*u_s(1,i+1)-u_s(1,i+2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta)) * ...
                ((-8*W_s(1,i)+9*W_s(1,i+1)-W_s(1,i+2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(1,i)/W_s(1,i))^2;
            
        elseif i==Ms
            de2_s(j,i) = ((u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((8*u_s(j,i)-9*u_s(j,i-1)+u_s(j,i-2))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2 + ...
                (u_s(j,i)/W_s(j,i))^2*((8*W_s(j,i)-9*W_s(j,i-1)+W_s(j,i-2))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta))^2 + ...
                u_s(j,i)/W_s(j,i)*((8*u_s(j,i)-9*u_s(j,i-1)+u_s(j,i-2))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta)) * ...
                ((8*W_s(j,i)-9*W_s(j,i-1)+W_s(j,i-2))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(j,i)/W_s(j,i))^2;
            
            de2_s(N,i) = ((3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((8*u_s(N,i)-9*u_s(N,i-1)+u_s(N,i-2))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2 + ...
                (u_s(N,i)/W_s(N,i))^2*((8*W_s(N,i)-9*W_s(N,i-1)+W_s(N,i-2))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta))^2 + ...
                u_s(N,i)/W_s(N,i)*((8*u_s(N,i)-9*u_s(N,i-1)+u_s(N,i-2))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta)) * ...
                ((8*W_s(N,i)-9*W_s(N,i-1)+W_s(N,i-2))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(N,i)/W_s(N,i))^2;
            
            de2_s(1,i) = ((-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((8*u_s(1,i)-9*u_s(1,i-1)+u_s(1,i-2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2 + ...
                (u_s(1,i)/W_s(1,i))^2*((8*W_s(1,i)-9*W_s(1,i-1)+W_s(1,i-2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta))^2 + ...
                u_s(1,i)/W_s(1,i)*((8*u_s(1,i)-9*u_s(1,i-1)+u_s(1,i-2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta)) * ...
                ((8*W_s(1,i)-9*W_s(1,i-1)+W_s(1,i-2))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(1,i)/W_s(1,i))^2;
            
        elseif i==2
            de2_s(j,i) = ((u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-4*u_s(j,i-1)+3*u_s(j,i)+u_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2 + ...
                (u_s(j,i)/W_s(j,i))^2*((-4*W_s(j,i-1)+3*W_s(j,i)+W_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta))^2 + ...
                u_s(j,i)/W_s(j,i)*((-4*u_s(j,i-1)+3*u_s(j,i)+u_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta)) * ...
                ((-4*W_s(j,i-1)+3*W_s(j,i)+W_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(j,i)/W_s(j,i))^2;
            
            de2_s(N,i) = ((3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-4*u_s(N,i-1)+3*u_s(N,i)+u_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2 + ...
                (u_s(N,i)/W_s(N,i))^2*((-4*W_s(N,i-1)+3*W_s(N,i)+W_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta))^2 + ...
                u_s(N,i)/W_s(N,i)*((-4*u_s(N,i-1)+3*u_s(N,i)+u_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta)) * ...
                ((-4*W_s(N,i-1)+3*W_s(N,i)+W_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(N,i)/W_s(N,i))^2;
            
            de2_s(1,i) = ((-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-4*u_s(1,i-1)+3*u_s(1,i)+u_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2 + ...
                (u_s(1,i)/W_s(1,i))^2*((-4*W_s(1,i-1)+3*W_s(1,i)+W_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta))^2 + ...
                u_s(1,i)/W_s(1,i)*((-4*u_s(1,i-1)+3*u_s(1,i)+u_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta)) * ...
                ((-4*W_s(1,i-1)+3*W_s(1,i)+W_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(1,i)/W_s(1,i))^2;
            
        elseif i==Ms-1
            de2_s(j,i) = ((u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-1*u_s(j,i-1)-3*u_s(j,i)+4*u_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2 + ...
                (u_s(j,i)/W_s(j,i))^2*((-1*W_s(j,i-1)-3*W_s(j,i)+4*W_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta))^2 + ...
                u_s(j,i)/W_s(j,i)*((-1*u_s(j,i-1)-3*u_s(j,i)+4*u_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta)) * ...
                ((-1*W_s(j,i-1)-3*W_s(j,i)+4*W_s(j,i+1))/(3*dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(j,i)/W_s(j,i))^2;
            
            de2_s(N,i) = ((3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-1*u_s(N,i-1)-3*u_s(N,i)+4*u_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2 + ...
                (u_s(N,i)/W_s(N,i))^2*((-1*W_s(N,i-1)-3*W_s(N,i)+4*W_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta))^2 + ...
                u_s(N,i)/W_s(N,i)*((-1*u_s(N,i-1)-3*u_s(N,i)+4*u_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta)) * ...
                ((-1*W_s(N,i-1)-3*W_s(N,i)+4*W_s(N,i+1))/(3*dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(N,i)/W_s(N,i))^2;
            
            de2_s(1,i) = ((-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((-1*u_s(1,i-1)-3*u_s(1,i)+4*u_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2 + ...
                (u_s(1,i)/W_s(1,i))^2*((-1*W_s(1,i-1)-3*W_s(1,i)+4*W_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta))^2 + ...
                u_s(1,i)/W_s(1,i)*((-1*u_s(1,i-1)-3*u_s(1,i)+4*u_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta)) * ...
                ((-1*W_s(1,i-1)-3*W_s(1,i)+4*W_s(1,i+1))/(3*dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(1,i)/W_s(1,i))^2;
            
        else
            de2_s(j,i) = ((u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((u(j,i)-u(j,i-1))/(dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))^2 + ...
                (u_s(j,i)/W_s(j,i))^2*((W(j,i)-W(j,i-1))/(dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta))^2 + ...
                u_s(j,i)/W_s(j,i)*((u(j,i)-u(j,i-1))/(dx) +...
                dzetadx_s(j,i)*(u_s(j+1,i)-u_s(j-1,i))/(2*dzeta))*...
                ((W(j,i)-W(j,i-1))/(dx) +...
                dzetadx_s(j,i)*(W_s(j+1,i)-W_s(j-1,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(j,i)/W_s(j,i))^2;
            
            de2_s(N,i) = ((3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((u(N,i)-u(N,i-1))/(dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))^2 + ...
                (u_s(N,i)/W_s(N,i))^2*((W(N,i)-W(N,i-1))/(dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta))^2 + ...
                u_s(N,i)/W_s(N,i)*((u(N,i)-u(N,i-1))/(dx) +...
                dzetadx_s(N,i)*(3*u_s(N,i)-4*u_s(N-1,i)+u_s(N-2,i))/(2*dzeta))*...
                ((W(N,i)-W(N,i-1))/(dx) +...
                dzetadx_s(N,i)*(3*W_s(N,i)-4*W_s(N-1,i)+W_s(N-2,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(N,i)/W_s(N,i))^2;
            
            de2_s(1,i) = ((-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2/(4*H_s(i)^2) + ...
                ((u(1,i)-u(1,i-1))/(dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))^2 + ...
                (u_s(1,i)/W_s(1,i))^2*((W(1,i)-W(1,i-1))/(dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta))^2 + ...
                u_s(1,i)/W_s(1,i)*((u(1,i)-u(1,i-1))/(dx) +...
                dzetadx_s(1,i)*(-3*u_s(1,i)+4*u_s(2,i)-u_s(3,i))/(2*dzeta))*...
                ((W(1,i)-W(1,i-1))/(dx) +...
                dzetadx_s(1,i)*(-3*W_s(1,i)+4*W_s(2,i)-W_s(3,i))/(2*dzeta)) + ...
                is_flowband/4*(u_s(1,i)/W_s(1,i))^2;
        end
        
    end
end
visc_s = 1/2*AGlen_s.^(-1/n).*(de2_s + de0).^((1-n)/(2*n)); % [Pa a]
visc = staggerX2main(visc_s);
strain_heat_s = 4*visc_s.*de2_s; % [Pa a-1]
strain_heat = staggerX2main(strain_heat_s); % [Pa a-1]
end