% practice 2
% Cartesian > Keplerian
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.209-215
clear;
% Initial Value
POSITION_KM  = [ -6045  -3490  2500];
VELOCITY_KMS = [ -3.457 6.618 2.533];

% Constant
MU_KM3S2  = 398600;

% Calculate
[h, i, OMEGA, e, omega, theta]  = cartesianToKeplerian(POSITION_KM, VELOCITY_KMS, MU_KM3S2);
fprintf("h: %.0f \ni: %.1f \nƒ¶: %.1f \ne: %.4f \nƒÖ: %.2f \nƒÆ: %.2f\n", h, i, OMEGA, e, omega, theta);

%%
%Cartesian > Keplerian
function [hNorm, incliDeg, OMEGADeg, eNorm, omegaDeg, thetaDeg]=cartesianToKeplerian(R,V,MU)
KUnitVector = [0 0 1];
% Calculation
distance = norm(R);
radialV  = dot(R, V) / distance;

% h: angular momentum
hVector  = cross(R, V);
hNorm    = norm(hVector);

% i: inclination
incliRad = acos(hVector(3) / hNorm);

% OMEGA: right ascension of the ascending node
NVector  = cross(KUnitVector, hVector);
OMEGARad = calcOMEGA(NVector);

% e: eccentricity
eVector  = 1 / MU * (cross(V, hVector) - MU * R / distance);
eNorm    = norm(eVector);

% omega: argument of perigee
omegaRad = calcomethe(NVector, eVector, eVector(3));

% theta: true anomaly
thetaRad = calcomethe(eVector, R, radialV);

% Rad > Deg
incliDeg = radToDeg(incliRad);
OMEGADeg = radToDeg(OMEGARad);
omegaDeg = radToDeg(omegaRad);
thetaDeg = radToDeg(thetaRad);
    %%
    % Radian >> Degree
    function deg = radToDeg(rad)
    deg = rad * 180 / pi;
    end
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