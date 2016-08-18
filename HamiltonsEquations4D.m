function out = hamiltonsEquations4D(time, p)
% Inputs:
% * time: a 1x2 row vector containing [initialTime, finalTime].
% * par : a 1x32 row vector containing position, momemtum, mass, and charge
%   parameters of the four particles. The particles' coordinate components are
%   given as [x1, y1, z1, ..., z4, px1, ..., pz4, m1, ..., m4, q1, ..., q4].
%   Had to do it this way because ode45 wants all the parameters in a single
%   vector.
% Output:
% * out : an nx25 array where each row contains
%         [time, position[1x12], momentum[1x12]] in the same format as par. In
%         practice, only the final row is utilized to evaluate the final
%         conditions of the system but you can use all the rows to study the
%         time evolution of the Coulomb explosion.
% Notes: All units are SI.

% Constants.
amu = 1.66053886e-27; % [kg], 1 atomic mass unit
e   = 1.60217646e-19; % [C], 1 elementary charge
k   = 8.987551e9;     % [N m^2 C^-2], electrostatic constant

% Masses and charges
m1 = amu*p(25); m2 = amu*p(26); m3 = amu*p(27); m4 = amu*p(28);
q1 = e*p(29);   q2 = e*p(30);   q3 = e*p(31);   q4 = e*p(32);

% Calculate the distance between ions. Note that this quantity does not preserve
% vector direction. All in [m].
r12 = sqrt((p(1)-p(4))^2  + (p(2)-p(5))^2  + (p(3)-p(6))^2);
r13 = sqrt((p(1)-p(7))^2  + (p(2)-p(8))^2  + (p(3)-p(9))^2);
r14 = sqrt((p(1)-p(10))^2 + (p(2)-p(11))^2 + (p(3)-p(12))^2);
r23 = sqrt((p(4)-p(7))^2  + (p(5)-p(8))^2  + (p(6)-p(9))^2);
r24 = sqrt((p(4)-p(10))^2 + (p(5)-p(11))^2 + (p(6)-p(12))^2);
r34 = sqrt((p(7)-p(10))^2 + (p(8)-p(11))^2 + (p(9)-p(12))^2);

% pDot is a column vector with components
% [vx1; vy1; vz1; ...; vz4; p'x1; ... p'z4].
% These quantities are produced by taking the first derivative of the
% Hamiltonian with respect to the appropriate variable.

v1 = p(13:15)' / m1;
v2 = p(16:18)' / m2;
v3 = p(19:21)' / m3;
v4 = p(22:24)' / m4;

% pDot in this case is just the pairwise Coulomb force so just call it F_ij.
F12 = (k*q1*q2 / r12^3) * [(p(1)-p(4));  (p(2)-p(5));  (p(3)-p(6));];
F13 = (k*q1*q3 / r13^3) * [(p(1)-p(7));  (p(2)-p(8));  (p(3)-p(9));];
F14 = (k*q1*q4 / r14^3) * [(p(1)-p(10)); (p(2)-p(11)); (p(3)-p(12));];
F23 = (k*q2*q3 / r23^3) * [(p(4)-p(7));  (p(5)-p(8));  (p(6)-p(9));];
F24 = (k*q2*q4 / r24^3) * [(p(4)-p(10)); (p(5)-p(11)); (p(6)-p(12));];
F34 = (k*q3*q4 / r34^3) * [(p(7)-p(10)); (p(8)-p(11)); (p(9)-p(12));];

% TODO: Double check the plus/minus signs here!
p1Dot =  F12 + F13 + F14;
p2Dot = -F12 + F23 + F24;
p3Dot = -F13 - F23 + F34;
p4Dot = -F14 - F24 - F34;

pDot = [v1; v2; v3; v4; p1Dot; p2Dot; p3Dot; p4Dot;];

% Output must be of the same form as the input parameters p so we append the
% mass and charge values that were initially there.
out = [parDot; p(25:32)];
end
