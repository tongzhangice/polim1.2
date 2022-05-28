function set_staggered_grid()

global M Ms N dx zeta xi xi_s
global hS hB H hS_s hB_s H_s
global dhSdx_s dhSdx dhBdx_s dhBdx dzetadx_s dzetadx
global beta2_s

xi_s = main2staggerX(xi);
hS_s = main2staggerX(hS);
hB_s = main2staggerX(hB);
H_s = main2staggerX(H);

dhSdx_s = zeros(1,Ms);
dhSdx_s(2:M) = (hS(2:end)-hS(1:end-1))/dx;
dhSdx_s(1) = (hS_s(2)-hS_s(1))/dx;
dhSdx_s(end) = (hS_s(end)-hS_s(end-1))/dx;
dhSdx = staggerX2main(dhSdx_s);
% dhSdx = [dhSdx_s(1), (dhSdx_s(2:end-2)+dhSdx_s(3:end-1))/2, dhSdx_s(end)];

dhBdx_s = zeros(1,Ms);
dhBdx_s(2:M) = (hB(2:end)-hB(1:end-1))/dx;
dhBdx_s(1) = (hB_s(2)-hB_s(1))/dx;
dhBdx_s(end) = (hB_s(end)-hB_s(end-1))/dx;
dhBdx = staggerX2main(dhBdx_s);
% dhBdx = [dhBdx_s(1), (dhBdx_s(2:end-2)+dhBdx_s(3:end-1))/2, dhBdx_s(end)];

dzetadx_s = -((1-zeta)*dhBdx_s + zeta*dhSdx_s)./(ones(N,1)*H_s);
dzetadx = staggerX2main(dzetadx_s);

% dzetadx = -((1-zeta)*dhBdx + zeta*dhSdx)./(ones(N,1)*H);
% dzetadx_s = main2staggerX(dzetadx);

% beta2 = 1e3 + 1e3*sin(2*pi/xi(end)*xi);
% beta2 = 1e4 + 1e4*sin(2*pi/xi(end)*xi);
% beta2 = 1e5*ones(1,M);
% beta2_s = main2staggerX(beta2);

end