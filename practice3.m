% practice 3
% Keplerian > Cartesian
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.231-233

% Orbital Elements
h_KM2S = 80000;
e = 1.4;
incli_DEG = 30;
OMEGA_DEG = 40;
omega_DEG = 60;
theta_DEG = 30;

% Constant
MU_KM3S2  = 398600;

StateVector = KeplerianToCartesian(h_KM2S, e, incli_DEG, OMEGA_DEG, omega_DEG, theta_DEG, MU_KM3S2);
fprintf(" %.0f\n  %.0f\n  %.0f\n%.2f\n%.3f\n %.3f\n",StateVector);

%%
function StateVector = KeplerianToCartesian(h_KM2S, e, incli_DEG, OMEGA_DEG, omega_DEG, theta_DEG, MU_KM3S2)
incli_RAD = DegToRad(incli_DEG);
OMEGA_RAD = DegToRad(OMEGA_DEG);
omega_RAD = DegToRad(omega_DEG);
theta_RAD = DegToRad(theta_DEG);

position_per = h_KM2S ^ 2 / MU_KM3S2 / (1 + e * cos(theta_RAD)) * [cos(theta_RAD); sin(theta_RAD); 0];
velocity_per = MU_KM3S2 / h_KM2S * [-sin(theta_RAD); e + cos(theta_RAD); 0];
Q_pergeo = PerToGeo(OMEGA_RAD, incli_RAD, omega_RAD);
position_geo = Q_pergeo * position_per;
velocity_geo = Q_pergeo * velocity_per;
StateVector = [position_geo, velocity_geo];
    %%
    function rad = DegToRad(deg)
    rad = deg / 180 * pi;
    end
    %%
    function Qpg = PerToGeo(OMEGA,incli,omega)
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