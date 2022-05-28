function [Hnp1] = get_evolution_continuity(Hn, u, SMB, dt_u)
% The continuity equation is converted into a diffusion equation (Pattyn 2002,
% Pimentel 2011).
% Discretization of the continuity equation uses the semi-implicit metHnd

% Inputs:
% u: horizontal velocity field
% SMB: mass balance [m a-1]
% dt: time step [a]

% Outputs:
% Hn: ice thickness

global M dzeta dx hS hB W dhSdx

type_formulation = 2;
Wsurf = W(end,:);

uav = dzeta*(sum(u(1:end-1,:))+sum(u(2:end,:)))/2; % [1 by M]
uav_s = main2staggerX(uav);

diffu = zeros(1,M);
for i = 2:M-1
    diffu(i) = -uav_s(i+1).*Hn(i)./(dhSdx(i)+1e-10);
end
diffu(1) = -uav_s(2).*Hn(1)./((hS(2)-hS(1))/dx+1e-10);
diffu(M) = -uav_s(M+1).*Hn(M)./((hS(M)-hS(M-1))/dx+1e-10);

diffu_mid = (diffu(1:end-1)+diffu(2:end))/2; % 1 by M-1
LT1 = zeros(M,1);
LT2 = zeros(M,1);
LT3 = zeros(M,1);
RT = zeros(M,1);

switch type_formulation
    case 1
        % D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]
        % forward difference for width variation
        for i = 2:M-1
            LT1(i) = -diffu_mid(i-1)/(dx^2);
            LT2(i) = 1/dt_u + (diffu_mid(i) + diffu_mid(i-1))/(dx^2) +...
                diffu_mid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);
            LT3(i) = -diffu_mid(i)/(dx^2) -...
                diffu_mid(i)*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2);
            
            RT(i) = Hn(i)/dt_u + diffu_mid(i)*(hB(i+1)-hB(i))/(dx^2) -...
                diffu_mid(i-1)*(hB(i)-hB(i-1))/(dx^2) +...
                diffu_mid(i)*(hB(i+1)-hB(i))*(Wsurf(i+1)-Wsurf(i))/Wsurf(i)/(dx^2) +...
                SMB(i);
        end
    case 2
        % D = -uav * H / [(s(i+1) - s(i-1))/(2*dx)]
        % forward difference for width variation
        for i = 2:(M-1)
            LT1(i) = -diffu_mid(i-1)/(dx^2) +...
                diffu_mid(i)*(Wsurf(i+1)-Wsurf(i-1))/(4*dx^2*Wsurf(i));
            LT2(i) = 1/dt_u + (diffu_mid(i) + diffu_mid(i-1))/(dx^2);
            LT3(i) = -diffu_mid(i)/(dx^2) -...
                diffu_mid(i)*(Wsurf(i+1)-Wsurf(i-1))/(4*dx^2*Wsurf(i));
            RT(i) = Hn(i)/dt_u + diffu_mid(i)*(hB(i+1)-hB(i))/(dx^2) -...
                diffu_mid(i-1)*(hB(i)-hB(i-1))/(dx^2) +...
                diffu_mid(i)*(Wsurf(i+1)-Wsurf(i-1))*(hB(i+1)-hB(i-1))/(4*dx^2*Wsurf(i)) +...
                SMB(i);
        end
end

% boundary condition at glacier head
LT2(1) = 1; RT(1) = -uav_s(2)*Hn(1)*dt_u/dx + SMB(1);

% boundary condition at glacier terminus
LT2(end) = 1; RT(end) = uav_s(end)*Hn(end)*dt_u/dx + SMB(end);

LT = spdiags([[LT1(2:end); 0], LT2, [0; LT3(1:end-1)]], [-1, 0, 1], M, M);
Hnp1 = LT\RT; % [M by 1]
Hnp1 = Hnp1'; % [1 by M]
end
