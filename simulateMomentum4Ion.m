% simulateMomenta takes in
% * geometries: a matrix where each row is of the form [r_12 r_23 theta].
%   r_12 and r_23 should be given in SI units [m] and theta in [deg].
% * masses:     a vector [m1 m2 m3] with the atomic masses in amu.
% * charges:    a vector [q1 q2 q3] with the atomic charges in units of e.
function out = simulateMomentum4Ion(ZMatrix, masses, charges)
  geometry = ZMatrixToCartesian(ZMatrix);
  p_0 = zeros(1,12);
  p = timeEvolveMomenta4Ion([g p_0], masses, charges);
  out = rotateMomentum4Ion(p);
end

% IonVelocities: Calculate velocities of ion given initial parameters
% Usage: IonVelocities(InitialConditions) where InitialConditions is a vector
% of the form [x1, y1, z1, x2, ..., z3].
function out = timeEvolveMomenta4Ion(initialConditions, masses, charges)
  options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
  [t,y] = ode45('hamiltonsEquations4Ion', [0 1e-11], [initialConditions masses charges], options);
  out = y(size(t,1), 13:24);
end
