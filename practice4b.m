%%
% practice 4b
% Initial value > State vector
% Reference:
% îºógèróY, "É~ÉbÉVÉáÉìâêÕÇ∆ãOìπê›åvÇÃäÓëb", pp.37-39.
% 4aÇ∆ç≈å„ÇÃä÷êîÇÃÇ›àŸÇ»ÇÈ
%%
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;
t0 = 0;
% Constant
MU_KM3S2 = 398600;

stateVector = method2(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2);
fprintf("  %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);
%%
function stateVector = method2(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2)
%Calculation
distance = norm(POSITION_KM);
energy   = dot(VELOCITY_KMS, VELOCITY_KMS) / 2 - MU_KM3S2 / distance;

hVector  = cross(POSITION_KM, VELOCITY_KMS);
hNorm    = norm(hVector);
PVector  = (-MU_KM3S2 / distance) * POSITION_KM - cross(hVector, VELOCITY_KMS);
PNorm    = norm(PVector);

semiMajorAxis   = - MU_KM3S2 / (2 * energy);
eccentricity    = PNorm / MU_KM3S2;
semiLatusRectum = (hNorm ^ 2) / MU_KM3S2;
if energy < 0
    ecosz0  = 1 - distance / semiMajorAxis;
    esinz0  = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiMajorAxis);
    z0      = atan(esinz0 / ecosz0);
    t_pi    = t0 - sqrt(semiMajorAxis ^ 3 / MU_KM3S2) * (z0 - esinz0);
elseif energy > 0
    esinhz0 = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * (-semiMajorAxis));
    z0      = asinh(dot(POSITION_KM, VELOCITY_KMS) / (eccentricity * sqrt(MU_KM3S2 * (-semiMajorAxis))));
    t_pi    = t0 - sqrt((-semiMajorAxis) ^ 3 / MU_KM3S2) * (esinhz0 - z0);
else
    z0      = dot(POSITION_KM, VELOCITY_KMS) / sqrt(MU_KM3S2 * semiLatusRectum);
    t_pi    = t0 - 0.5 * sqrt(semiLatusRectum ^ 3 / MU_KM3S2) * (z0 + z0 ^ 3 / 3);
end

z = calcKepler(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, z0, DELTA_TIME_S, t_pi);
stateVector = LagtoState(MU_KM3S2, semiMajorAxis, semiLatusRectum, energy, POSITION_KM, VELOCITY_KMS,...
                         z, z0, DELTA_TIME_S, t0);
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
    function stateVector = LagtoState(mu, a, p, energy, R0, V0, z, z0, t, t0)
    r0 = norm(R0);
    delta_z = z - z0;
    if energy < 0
        f = 1 - a / r0 * (1 - cos(delta_z));
        g = t - t0 - sqrt(a ^ 3 / mu) * (delta_z - sin(delta_z));
    elseif energy > 0
        f = 1 - a / r0 * (1 - cosh(delta_z));
        g = t - t0 + sqrt((- a) ^ 3 / mu) * (delta_z - sinh(delta_z));
    else
        f = 1 - p / (2 * r0) * (delta_z ^ 2);
        g = t - t0 - 1 / 6 * sqrt(p ^ 3 / mu) * (delta_z ^ 3);
    end
    R = f * R0 + g * V0;
    r = norm(R);
    if energy < 0
        fd = - sqrt(mu * a) / (r * r0) * sin(delta_z);
        gd = 1 - a / r * (1 - cos(delta_z));
    elseif energy > 0
        fd = - sqrt(mu * (- a)) / (r * r0) * sinh(delta_z);
        gd = 1 - a / r * (1 - cosh(delta_z));
    else
        fd = - sqrt(mu * p) / (r * r0) * delta_z;
        gd = 1 - p / (2 * r) * (delta_z ^ 2);
    end
    V = fd * R0 + gd * V0;
    stateVector = [R, V];
    end
end