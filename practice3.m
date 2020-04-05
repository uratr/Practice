% practice 3
% Keplerian > Cartesian
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.231-233
clear;
% Orbital Elements
h_KM2S = 80000;
e = 1.4;
incli_DEG = 30;
OMEGA_DEG = 40;
omega_DEG = 60;
theta_DEG = 30;

% Constant
MU_KM3S2  = 398600;

stateVector = keplerianToCartesian(h_KM2S, e, incli_DEG, OMEGA_DEG, omega_DEG, theta_DEG, MU_KM3S2);
fprintf(" %.0f\n  %.0f\n  %.0f\n%.2f\n%.3f\n %.3f\n",stateVector);

%%
function stateVector = keplerianToCartesian(h_KM2S, e, incli_DEG, OMEGA_DEG, omega_DEG, theta_DEG, MU_KM3S2)
incliRad = degToRad(incli_DEG);
OMEGARad = degToRad(OMEGA_DEG);
omegaRad = degToRad(omega_DEG);
thetaRad = degToRad(theta_DEG);

positionPer = h_KM2S ^ 2 / MU_KM3S2 / (1 + e * cos(thetaRad)) * [cos(thetaRad); sin(thetaRad); 0];
velocityPer = MU_KM3S2 / h_KM2S * [-sin(thetaRad); e + cos(thetaRad); 0];
Q_PerGeo    = perToGeo(OMEGARad, incliRad, omegaRad);
positionGeo = Q_PerGeo * positionPer;
velocityGeo = Q_PerGeo * velocityPer;
stateVector = [positionGeo, velocityGeo];
    %%
    function rad = degToRad(deg)
    rad = deg / 180 * pi;
    end
    %%
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
end