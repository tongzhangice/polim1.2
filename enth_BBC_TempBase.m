function [LT, RT] = enth_BBC_TempBase(LT1, LT2, LT3, RT, Epmp_i)
global N

LT1(1) = 0;
LT2(1) = 1;
LT3(1) = 0;
RT(1) = Epmp_i(1);

% RT(1) = Epmp_i(1) + 3.34e5*0.01;

LT = spdiags([[LT1(2:end);0], LT2, [0;LT3(1:end-1)]], [-1,0,1], N, N);

end