%%
% practice 6a
%%
% practice 4a
% Initial value > State vector
% Reference:
% 半揚俊雄, "ミッション解析と軌道設計の基礎", pp.37-39.
% 地球の扁平率を考慮していない．
% 時間変化が10sec程度であれば4acと同じ結果だが本条件では誤差が生じる．
%%
% 共通
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;
t0 = 0;
% Constant
MU_KM3S2 = 398600;
%{
%%
% Kepler Eq.
stateVector = method1(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2);
fprintf("Kepler Eq.\n %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);

%%
% RK4
h = 0.1;
t = t0;
V = VELOCITY_KMS;
X = POSITION_KM;
r = norm(POSITION_KM);
while t <= DELTA_TIME_S
     t  = t + 1;
     V_K1 = RK4_VEL(MU_KM3S2, r, X);
     V_K2 = RK4_VEL(MU_KM3S2, r, X + h / 2 * V_K1);
     V_K3 = RK4_VEL(MU_KM3S2, r, X + h / 2 * V_K2);
     V_K4 = RK4_VEL(MU_KM3S2, r, X + h * V_K3);
     V  = V + h / 6 * (V_K1 + 2 * V_K2 + 2 * V_K3 + V_K4);
     X_K1 = RK4_POS(V);
     X_K2 = RK4_POS(V + h / 2 * X_K1);
     X_K3 = RK4_POS(V + h / 2 * X_K2);
     X_K4 = RK4_POS(V + h * X_K3);
     X  = X + h / 6 * (X_K1 + 2 * X_K2 + 2 * X_K3 + X_K4);
end
fprintf("RK4:\n %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",X,V);
%}
%%
% ODE
tspan = [t0 DELTA_TIME_S];
X0 = POSITION_KM;
V0 = VELOCITY_KMS;

[t,y] = ode45(@myODE,tspan,[X0; V0]);
fprintf("ODE:\n %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",y(end,:));
function dy = myODE(t,Y)
mu = 398600;
POSITION_KM  = [ -3670 -3870 4400];
r = norm(POSITION_KM);
dy = [Y(4:6); - mu / (r ^ 3) * Y(1:3)];
end
%%
%{
% For Kepler Eq.
function stateVector = method1(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2)
%Calculation
distance = norm(POSITION_KM);
energy   = dot(VELOCITY_KMS, VELOCITY_KMS) / 2 - MU_KM3S2 / distance;

hVector  = cross(POSITION_KM, VELOCITY_KMS);
hNorm    = norm(hVector);
PVector  = (-MU_KM3S2 / distance) * POSITION_KM - cross(hVector, VELOCITY_KMS);
PNorm    = norm(PVector);

WUnitVector = hVector / hNorm;
PUnitVector = PVector / PNorm;
QUnitVector = cross(WUnitVector, PUnitVector);

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
stateVector = UnittoState(MU_KM3S2, semiMajorAxis, eccentricity, semiLatusRectum, energy, PUnitVector, QUnitVector, z);
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
    function stateVector = UnittoState(mu, a, e, p, energy, P, Q, z)
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
end
%%
% For RK4 velocity
function result = RK4_VEL(mu, r, x)
result = - mu / (r ^ 3) * x;
end
%%
% For RK4 position
function result = RK4_POS(v)
result = v;
end
%}