% Kadai 1
% Euler angle > direction cosine matrix

% Euler angle[deg]
% OMEGA = 
EULER_ANGLE_DEG = [40 30 60];
EULER_ANGLE_RAD = DegToRad(EULER_ANGLE_DEG);
% Direction cosine matrix
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