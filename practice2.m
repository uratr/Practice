% practice 2
% Cartesian > Keplerian
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.209-215

% Initial Value (Curtis)
POSITION_KM  = [ -6045  -3490  2500];
VELOCITY_KMS = [ -3.457 6.618 2.533];

% Constant
MU_KM3S2  = 398600;

% Calculate
[h, i, OMEGA, e, omega, theta]  = CartesianToKeplerian(POSITION_KM, VELOCITY_KMS, MU_KM3S2);
fprintf("h: %.0f \ni: %.1f \nƒ¶: %.1f \ne: %.4f \nƒÖ: %.2f \nƒÆ: %.2f\n", h, i, OMEGA, e, omega, theta);

%%
%Cartesian > Keplerian
function [hnorm, incli_deg, OMEGA_deg, enorm, omega_deg, theta_deg]=CartesianToKeplerian(R,V,MU)
K         = [0 0 1];
% Calculation
distance  = norm(R);
RadialV   = dot(R, V) / distance;

% h: angular momentum
hvector   = cross(R, V);
hnorm     = norm(hvector);

% i: inclination
incli_rad = acos(hvector(3) / hnorm);

% OMEGA: right ascension of the ascending node
Nvector   = cross(K, hvector);
OMEGA_rad = calc_OMEGA(Nvector);

% e: eccentricity
evector   = 1 / MU * (cross(V, hvector) - MU * R / distance);
enorm     = norm(evector);

% omega: argument of perigee
omega_rad = calc_omethe(Nvector, evector, evector(3));

% theta: true anomaly
theta_rad = calc_omethe(evector, R, RadialV);

% Rad > Deg
incli_deg = RadToDeg(incli_rad);
OMEGA_deg = RadToDeg(OMEGA_rad);
omega_deg = RadToDeg(omega_rad);
theta_deg = RadToDeg(theta_rad);
%%
% Radian >> Degree
function deg = RadToDeg(rad)
deg = rad * 180 / pi;
end
%%
% Calculate OMEGA
function OMEGA = calc_OMEGA(N)
if N(2) >= 0
    OMEGA = acos(N(1) / norm(N));
else
    OMEGA = 2 * pi - acos(N(1) / norm(N));
end
end
%%
% calculate omega & theta
function angle = calc_omethe(A, B, c)
if c >= 0
    angle = acos(dot(A / norm(A), B / norm(B)));
else
    angle = 2 * pi - acos(dot(A / norm(A), B / norm(B)));
end
end
end