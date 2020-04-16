%%
% practice 4c
% Initial value > State vector
% Reference:
% ”¼—gr—Y, "ƒ~ƒbƒVƒ‡ƒ“‰ðÍ‚Æ‹O“¹ÝŒv‚ÌŠî‘b", pp.74-76.
% universal variable
%%
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;
t0 = 0;
% Constant
MU_KM3S2 = 398600;

stateVector = method3(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2);
fprintf("  %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);
%%
function stateVector = method3(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, MU_KM3S2)
%Calculation
distance = norm(POSITION_KM);
alpha0   = 2 / distance - dot(VELOCITY_KMS, VELOCITY_KMS) / MU_KM3S2;
if alpha0 > 0
    x0 = alpha0 * sqrt(MU_KM3S2) * (DELTA_TIME_S - t0);
elseif alpha0 < 0
    x0 = 1 / sqrt(-alpha0) * log((2 * MU_KM3S2 *(sqrt(-alpha0)) ^ 3 * (DELTA_TIME_S - t0)) / ...
        (sqrt(MU_KM3S2) * (1 - alpha0 * distance) + sqrt(-alpha0) * dot(POSITION_KM, VELOCITY_KMS)));
else
    fprintf("ERROR\n");
end
x = calc_x(DELTA_TIME_S,t0,POSITION_KM,VELOCITY_KMS,alpha0,x0,MU_KM3S2);
stateVector = LagtoState(MU_KM3S2, POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, t0, x, alpha0);
    %%
    % calculate position & velocity
    function stateVector = LagtoState(mu, R0, V0, t, t0, x, alpha0)
    r0 = norm(R0);
    f = 1 - (x ^ 2) / r0 * funC(alpha0 * (x ^ 2));
    g = t - t0 - (x ^ 3) / sqrt(mu) * funS(alpha0 * (x ^ 2));
    R = f * R0 + g * V0;
    r = norm(R);
    fd = sqrt(mu) / (r * r0) * (alpha0 * (x ^ 3) * funS(alpha0 * (x ^ 2)) - x);
    gd = 1 - (x ^ 2) / r * funC(alpha0 * (x ^ 2));
    V = fd * R0 + gd * V0;
    stateVector = [R, V];
    end
%%
    function result = calc_x(t,t0,R0,V0,alpha0,x0,mu)
        r0   = norm(R0);
        xp    = x0;
        for ix = 1:100
            z    = alpha0 * (xp ^ 2);
            S_z  = funS(z);
            C_z  = funC(z);
            F_x  = (1 - alpha0 * r0) * (xp ^ 3) * S_z +...
                   dot(R0, V0) / sqrt(mu) * (xp ^ 2) * C_z +...
                   r0 * xp - sqrt(mu) * (t - t0);
            dFdx = (1 - alpha0 * r0) * (xp ^ 2) * C_z +...
                   dot(R0, V0) / sqrt(mu) * xp * (1 - z * S_z) + r0;
            xpp  = xp - F_x / dFdx;
            if (abs (xpp - xp) > 1e-8) && (ix < 100)
                xp = xpp;
            elseif (abs (xpp - xp) > 1e-8) && (ix >= 100)
                fprintf("ERROR_x\n");
            else
                break
            end
        end
        result = xp;
    end
%%
    function result = funS(x)
        if x > 0
            result = (sqrt(x) - sin(sqrt(x))) / ((sqrt(x)) ^ 3);
        elseif x < 0
            result = (sinh(sqrt(-x)) - sqrt(-x)) / ((sqrt(-x)) ^ 3);
        else
            result  = 1 / factorial(3);
        end
    end
%%
    function result = funC(x)
        if x > 0
            result = (1 - cos(sqrt(x))) / x;
        elseif x < 0
            result = (cosh(sqrt(-x)) - 1) / (-x);
        else
            result  = 1 / factorial(2);
        end
    end
end