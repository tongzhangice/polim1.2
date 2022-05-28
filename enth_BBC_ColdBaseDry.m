function [LT, RT] = enth_BBC_ColdBaseDry(LT1, LT2, LT3, RT, H_icol, para)
global N dzeta

Qgeo = para.Qgeo;
Cp = para.Cp;
kc = para.kc;

LT1(1) = 0;
LT2(1) = -1/dzeta;
LT3(1) = 1/dzeta;
RT(1) = -Qgeo*H_icol*Cp/kc;
LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);

end