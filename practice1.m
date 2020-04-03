% practice 1
% Euler angle > direction cosine matrix(DCM)
% Reference: Curtis. H, "Orbital Mechanics for engineering students. 2nd ed". pp.230-232

% Euler angle[deg]
% [OMEGA i omega]
% OMEGA: the right ascension of the ascending node
% incli: the orbital inclination angle
% omega: the argument of perigee
EULER_ANGLE_DEG = [40 30 60];
EULER_ANGLE_RAD = DegToRad(EULER_ANGLE_DEG);

% DCM
Q = DirectionCosMat(EULER_ANGLE_RAD(1), EULER_ANGLE_RAD(2), EULER_ANGLE_RAD(3));
disp(Q)

%%
function Rad = DegToRad(Deg)
Rad = Deg / 180 * pi;
end
%%
function Qmat = DirectionCosMat(OMEGA,incli,omega)
Qmat = [ -sin(OMEGA) * cos(incli) * sin(omega) + cos(OMEGA) * cos(omega)...
          cos(OMEGA) * cos(incli) * sin(omega) + sin(OMEGA) * cos(omega)...
          sin(incli) * sin(omega);
         -sin(OMEGA) * cos(incli) * cos(omega) - cos(OMEGA) * sin(omega)...
          cos(OMEGA) * cos(incli) * cos(omega) - sin(OMEGA) * sin(omega)...
          sin(incli) * cos(omega);
          sin(OMEGA) * sin(incli)...
         -cos(OMEGA) * sin(incli)...
          cos(incli)];
end