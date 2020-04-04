% practice 1
% Euler angle > direction cosine matrix(DCM)
% Reference: 
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.230-232

% Euler angle[deg]
% [OMEGA incli omega]
% OMEGA: the right ascension of the ascending node
% incli: the orbital inclination angle
% omega: the argument of perigee
EULER_ANGLE_DEG = [40 30 60];
eulerAngleRad = degToRad(EULER_ANGLE_DEG);

% DCM
Q = geoToPer(eulerAngleRad(1), eulerAngleRad(2), eulerAngleRad(3));
disp(Q)

%%
function rad = degToRad(deg)
rad = deg / 180 * pi;
end
%%
function Qgp = geoToPer(OMEGA,incli,omega)
Qgp = [ -sin(OMEGA) * cos(incli) * sin(omega) + cos(OMEGA) * cos(omega)...
         cos(OMEGA) * cos(incli) * sin(omega) + sin(OMEGA) * cos(omega)...
         sin(incli) * sin(omega);
        -sin(OMEGA) * cos(incli) * cos(omega) - cos(OMEGA) * sin(omega)...
         cos(OMEGA) * cos(incli) * cos(omega) - sin(OMEGA) * sin(omega)...
         sin(incli) * cos(omega);
         sin(OMEGA) * sin(incli)...
        -cos(OMEGA) * sin(incli)...
         cos(incli)];
end