% practice 4b
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

hVector  = cross(POSITION_KM, VELOCITY_KMS);
hNorm    = norm(hVector);
PVector  = - MU_KM3S2 / distance * POSITION_KM - cross(hVector, VELOCITY_KMS);
PNorm    = norm(PVector);

WUnitVector = hVector / hNorm;
PUnitVector = PVector / PNorm;
QUnitVector = cross(WUnitVector, PUnitVector);

semiMajorAxis   = - MU_KM3S2 / (2 * energy);
eccentricity    = PNorm / MU_KM3S2;
semiLatusRectum = hNorm ^ 2 / MU_KM3S2;

period = 2 * pi / sqrt(MU_KM3S2) * semiMajorAxis ^ 1.5;
meanMotion = 2 * pi / period;
if energy < 0
    ecosz0 = 1 - distance / semiMajorAxis;
    esinz0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiMajorAxis);
    z0 = atan(esinz0 / ecosz0);
    t0 = (z0 - eccentricity * sin(z0)) / meanMotion;
    t_pi = t0 - sqrt(semiMajorAxis ^ 3 / MU_KM3S2) * (z0 - esinz0);
elseif energy > 0
    esinhz0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * (-semiMajorAxis));
    z0 = asinh(dot(POSITION_KM, VELOCITY_KMS) / (eccentricity * sqrt(MU_KM3S2 * (-semiMajorAxis))));
    t_pi = t0 - sqrt((-semiMajorAxis) ^ 3 / MU_KM3S2) * (esinhz0 - z0);
else
    z0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiLatusRectum);
    t_pi = t0 - 0.5 * sqrt(semiLatusRectum ^ 3 / MU_KM3S2) * (z0 + z0 ^ 3 / 3);
end
z = calcKepler(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, energy, z0, DELTA_TIME_S, t_pi);
stateVector = KeptoState(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, energy, PUnitVector, QUnitVector, z);
fprintf("  %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);
%%
% Kepler Eq
function zipp = calcKepler(mu, a, e, p, energy, z0, t, t_pi)
nKepler = 100;  % number of iteration of Kepler equation
TOLERANCE = 1e-8;
zipp = z0;
for iKepler = 1:nKepler
    z = zipp;
    if energy < 0
        fz    = z - e * sin(z) - sqrt(mu / a ^ 3) * (t - t_pi);
        fzdot = 1 - e * cos(z);
    elseif energy > 0
        fz    = e * sinh(z) - z - sqrt(mu / (- a) ^ 3) * (t - t_pi);
        fzdot = e * cosh(z) - 1;
    else
        fz    = z + z ^ 3 / 3 - 2 * sqrt(mu / p ^ 3) * (t - t_pi);
        fzdot = 1 + z ^ 2;
    end
    zipp = z - fz / fzdot;
    if abs(zipp) > TOLERANCE
        continue
    else
        break
    end
end
end
%%
% calculate position & velocity
function stateVector = KeptoState(mu, a, e, p, energy, P, Q, z)
if energy < 0
    R = a * (cos(z) - e) * P + sqrt(a * p) * sin(z) * Q;
    r = norm(R);
    V = - sqrt(mu * a) / r * sin(z) * P + sqrt(mu * p) / r * cos(z) * Q;
elseif energy > 0
    R = - a * (e - cosh(z)) * P + sqrt(- a * p) * sinh(z) * Q;
    r = norm(R);
    V = - sqrt(mu * (- a)) / r * sinh(z) * P + sqrt(mu * p) / r * cosh(z) * Q;
else
    R = p / 2 * (1 - z ^ 2) * P + p * z * Q;
    r = norm(R);
    V = - sqrt(mu * p) / r * z * Q + sqrt(mu * p) / r * Q;
end
stateVector = [R, V];
end