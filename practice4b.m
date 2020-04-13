%%
% practice 4b
% Initial value > State vector 
% Use Lagrange's coefficient
% Reference:
% Curtis. H, "Orbital Mechanics for engineering students. 2nd ed", pp.233-244
% consider the earthÅfs oblateness
%%
clear;
% Initial Value
POSITION_KM  = [ -3670 -3870 4400];
VELOCITY_KMS = [   4.7  -7.4    1];
DELTA_TIME_S = 96 * 3600;

distance = norm(POSITION_KM);
spead = norm(VELOCITY_KMS);
