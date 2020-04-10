% practice 4b
% Initial value > State vector
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.123
% Use Lagrange's coefficient
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;
t0 = 0;
% Constant
MU_KM3S2 = 398600;
elems = cartesianToKeplerian(POSITION_KM,VELOCITY_KMS,MU_KM3S2);
%{
PUnitVector = [ -sin(elems(3)) * cos(elems(2)) * sin(elems(5)) + cos(elems(3)) * cos(elems(5));
                 cos(elems(3)) * cos(elems(2)) * sin(elems(5)) + sin(elems(3)) * cos(elems(5));
                 sin(elems(2)) * sin(elems(5))];
QUnitVector = [ -sin(elems(3)) * cos(elems(2)) * cos(elems(5)) - cos(elems(3)) * sin(elems(5));
                 cos(elems(3)) * cos(elems(2)) * cos(elems(5)) - sin(elems(3)) * sin(elems(5));
                 sin(elems(2)) * cos(elems(5))];
%}
PUnitVector = [ -sin(elems(3)) * cos(elems(2)) * sin(elems(5)) + cos(elems(3)) * cos(elems(5)) ...
                -sin(elems(3)) * cos(elems(2)) * cos(elems(5)) - cos(elems(3)) * sin(elems(5)) ...
                 sin(elems(2)) * cos(elems(5))];
QUnitVector = [  cos(elems(3)) * cos(elems(2)) * sin(elems(5)) + sin(elems(3)) * cos(elems(5)) ...
                 cos(elems(3)) * cos(elems(2)) * cos(elems(5)) - sin(elems(3)) * sin(elems(5)) ...
                 sin(elems(3)) * sin(elems(2))];             
             
distance = norm(POSITION_KM);
energy   = dot(VELOCITY_KMS, VELOCITY_KMS) / 2 - MU_KM3S2 / distance;

hVector  = cross(POSITION_KM, VELOCITY_KMS);
hNorm    = norm(hVector);
PVector  = (-MU_KM3S2 / distance) * POSITION_KM - cross(hVector, VELOCITY_KMS);
PNorm    = norm(PVector);

WUnitVector = hVector / hNorm;
PUnitVector1 = PVector / PNorm;
QUnitVector1 = cross(WUnitVector, PUnitVector);

semiMajorAxis   = - MU_KM3S2 / (2 * energy);
eccentricity    = PNorm / MU_KM3S2;
semiLatusRectum = (hNorm ^ 2) / MU_KM3S2;

period = 2 * pi / sqrt(MU_KM3S2) * semiMajorAxis ^ 1.5;
meanMotion = 2 * pi / period;
if energy < 0
    ecosz0 = 1 - distance / semiMajorAxis;
    esinz0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiMajorAxis);
    z0 = atan(esinz0 / ecosz0);
    t_pi = t0 - sqrt(semiMajorAxis ^ 3 / MU_KM3S2) * (z0 - esinz0);
elseif energy > 0
    esinhz0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * (-semiMajorAxis));
    z0 = asinh(dot(POSITION_KM, VELOCITY_KMS) / (eccentricity * sqrt(MU_KM3S2 * (-semiMajorAxis))));
    t_pi = t0 - sqrt((-semiMajorAxis) ^ 3 / MU_KM3S2) * (esinhz0 - z0);
else
    z0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiLatusRectum);
    t_pi = t0 - 0.5 * sqrt(semiLatusRectum ^ 3 / MU_KM3S2) * (z0 + z0 ^ 3 / 3);
end
z = calcKepler(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, z0, DELTA_TIME_S, t_pi);

fprintf("%d,   %d,   %d\n",z,semiMajorAxis,eccentricity);
stateVector = KeptoState(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, energy, PUnitVector, QUnitVector, z);
fprintf("  %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);
%%
%Cartesian > Keplerian
function orbitalElementsRad=cartesianToKeplerian(R, V, MU)
K         = [0 0 1];
% Calculation
distance  = norm(R);
radialV   = dot(R, V) / distance;

% h: angular momentum
hVector   = cross(R, V);
hNorm     = norm(hVector);

% i: inclination
incliRad = acos(hVector(3) / hNorm);

% OMEGA: right ascension of the ascending node
NVector   = cross(K, hVector);
OMEGARad = calcOMEGA(NVector);

% e: eccentricity
eVector   = 1 / MU * (cross(V, hVector) - MU * R / distance);
eNorm     = norm(eVector);

% omega: argument of perigee
omegaRad = calcomethe(NVector, eVector, eVector(3));

% theta: true anomaly
thetaRad = calcomethe(eVector, R, radialV);
orbitalElementsRad = [hNorm, incliRad, OMEGARad, eNorm, omegaRad, thetaRad];
    %%
    % Calculate OMEGA
    function OMEGA = calcOMEGA(N)
    if N(2) >= 0
        OMEGA = acos(N(1) / norm(N));
    else
        OMEGA = 2 * pi - acos(N(1) / norm(N));
    end
    end
    %%
    % calculate omega & theta
    function angle = calcomethe(A, B, c)
    if c >= 0
        angle = acos(dot(A / norm(A), B / norm(B)));
    else
        angle = 2 * pi - acos(dot(A / norm(A), B / norm(B)));
    end
    end
end
%%
% Kepler Eq
function zipp = calcKepler(mu, a, e, p, z0, t, t_pi)
nKepler = 100;  % number of iteration of Kepler equation
TOLERANCE = 1e-8;
zipp = z0;
for iKepler = 1:nKepler
    zi = zipp;
    if e < 1
        fz    = zi - e * sin(zi) - sqrt(mu / (a ^ 3)) * (t - t_pi);
        fzdot = 1 - e * cos(zi);
    elseif e > 1
        fz    = e * sinh(zi) - zi - sqrt(mu / (- a) ^ 3) * (t - t_pi);
        fzdot = e * cosh(zi) - 1;
    else
        fz    = zi + zi ^ 3 / 3 - 2 * sqrt(mu / p ^ 3) * (t - t_pi);
        fzdot = 1 + zi ^ 2;
    end
    ratioi = fz / fzdot;
    if abs(ratioi) > TOLERANCE
        zipp = zi - ratioi;
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