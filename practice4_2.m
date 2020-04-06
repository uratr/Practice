% practice 4_1
% Initial value > State vector
% Reference:
% îºógèróY, "É~ÉbÉVÉáÉìâêÕÇ∆ãOìπê›åvÇÃäÓëb", pp.37-
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;
% Constant
MU_KM3S2 = 398600;

distance = norm(POSITION_KM);
energy   = dot(VELOCITY_KMS, VELOCITY_KMS) / 2 - MU_KM3S2 / distance;
hVector = cross(POSITION_KM, VELOCITY_KMS);
hNorm = norm(hVector);
WVector = hVector / hNorm;
semilatusRectum = hNorm ^ 2 / MU_KM3S2;

PVector = - MU_KM3S2 / distance * POSITION_KM - cross(hVector, VELOCITY_KMS); 