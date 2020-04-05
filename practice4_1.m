% practice 4_1
% Initial value > State vector (only ellipse)
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.240-244
% Use practice 2
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;

% Constant
MU_KM3S2 = 398600;
RADIUS_KM = 6378;  % radius of the earth

% Second Zonal Harmonics (Earth)
J2 = 1.08263e-3;

stateVector = ivtoState(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, MU_KM3S2, RADIUS_KM, J2);
fprintf("  %.0f\n  %.0f\n %.0f\n%.3f\n %.3f\n%.4f\n",stateVector);
%% 
% Initial value > State vector
function stateVector = ivtoState(POSITION_KM, VELOCITY_KMS, DELTA_TIME_S, MU_KM3S2, RADIUS_KM, J2)
elements = cartesianToKeplerian(POSITION_KM, VELOCITY_KMS, MU_KM3S2);
% ++++ %
% 1: h %
% 2: i %
% 3: ƒ¶ %
% 4: e %
% 5: ƒÖ %
% 6: ƒÆ %
% ++++ %

semimajorAxis = elements(1) ^ 2 / MU_KM3S2 / (1 - elements(4) ^ 2);
period = 2 * pi / sqrt(MU_KM3S2) * semimajorAxis ^ 1.5;
meanMotion = 2 * pi / period;
E0 = 2 * atan(sqrt((1 - elements(4)) / (1 + elements(4))) * tan(elements(6) / 2));
t0 = (E0 - elements(4) * sin(E0)) / meanMotion;
tf = t0 + DELTA_TIME_S;
np = tf / period;
tlast = tf - period * fix(np);
meanAnomaly = meanMotion * tlast;
Ei = calcKepler(E0, elements(4), meanAnomaly);
trueAnomaly =  2 * atan(sqrt((1 + elements(4)) / (1 - elements(4))) * tan(Ei / 2));
positionPer = elements(1) ^ 2 / (MU_KM3S2 * (1 + elements(4) * cos(trueAnomaly))) *...
              [ cos(trueAnomaly);
                sin(trueAnomaly);
                0];
velocityPer = [-MU_KM3S2 / elements(1) * sin(trueAnomaly);
                MU_KM3S2 / elements(1) * (elements(4) + cos(trueAnomaly));
                0];
OMEGAdot = - ((3 * sqrt(MU_KM3S2) * J2 * RADIUS_KM ^ 2) /...
              (2 * (1 - elements(4) ^ 2) ^ 2 * semimajorAxis ^ 3.5)) * cos(elements(2));
OMEGAlast = elements(3) + OMEGAdot * DELTA_TIME_S;
omegadot = - ((3 * sqrt(MU_KM3S2) * J2 * RADIUS_KM ^ 2) /...
              (2 * (1 - elements(4) ^ 2) ^ 2 * semimajorAxis ^ 3.5)) *...
              (2.5 * sin(elements(2)) ^ 2 - 2);
omegalast = elements(5) + omegadot * DELTA_TIME_S;
Qpg = perToGeo(OMEGAlast, elements(2), omegalast);
positionGeo = Qpg * positionPer;
velocityGeo = Qpg * velocityPer;
stateVector = [positionGeo, velocityGeo];
    %%
    % Kepler Eq
    function Ei = calcKepler(E0, e, meanAnomaly)
    nKepler = 100;  % number of iteration of Kepler equation
    TOLERANCE = 1e-8;
    Eipp = E0;
    for iKepler = 1:nKepler
        Ei = Eipp;
        fE = Ei - e * sin(Ei) - meanAnomaly;
        fEdot = 1 - e * cos(Ei);
        ratioi = fE / fEdot;
        if abs(ratioi) > TOLERANCE
            Eipp = Ei - ratioi;
        else
            break
        end
    end
    end
    %%
    % perifocal to geocentric
    function Qpg = perToGeo(OMEGA,incli,omega)
    Qgp = [ -sin(OMEGA) * cos(incli) * sin(omega) + cos(OMEGA) * cos(omega)...
             cos(OMEGA) * cos(incli) * sin(omega) + sin(OMEGA) * cos(omega)...
             sin(incli) * sin(omega);
            -sin(OMEGA) * cos(incli) * cos(omega) - cos(OMEGA) * sin(omega)...
             cos(OMEGA) * cos(incli) * cos(omega) - sin(OMEGA) * sin(omega)...
             sin(incli) * cos(omega);
             sin(OMEGA) * sin(incli)...
            -cos(OMEGA) * sin(incli)...
             cos(incli)];
    Qpg =  transpose(Qgp);
    end
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
end