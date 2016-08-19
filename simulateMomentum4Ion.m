% simulateMomenta takes in
% * geometries: a matrix where each row is of the form [r_12 r_23 theta].
%   r_12 and r_23 should be given in SI units [m] and theta in [deg].
% * masses:     a vector [m1 m2 m3] with the atomic masses in amu.
% * charges:    a vector [q1 q2 q3] with the atomic charges in units of e.
function out = simulateMomentum4Ion(ZMatrix, masses, charges)
  geometry = ZMatrixToCartesian(ZMatrix);
  p_0 = zeros(1,9);
  p = ionVelocities([g p_0], masses, charges);
  % out = [r_12 r_23 theta p(1:3) p(4:6) p(7:9)];

  % Extract each atom's momentum into 2D vectors in preparation to rotate.
  p_1 = p(1:3);
  p_2 = p(4:6);
  p_3 = p(7:9);
  p_4 = p(10:12);

  % Put each momentum into column vector form so we can use the usual matrix
  % multiplication.
  p_1 = p_1'; p_2 = p_2'; p_3 = p_3'; p_4 = p_4';

  % Rotate p_1 to point along the x-axis.
  % http://math.stackexchange.com/questions/114512/how-to-find-the-orthonormal-transformation-that-will-rotate-a-vector-to-the-x-ax
  p1x = p_1(1); p1y = p_1(2); p1z = p_1(3);
  theta = atan2(p1z, p1y);
  alpha = atan2(p1y*cos(theta) - p1z*sin(theta), p1x);
  R = Rz(-alpha)*Rx(-theta);
  p_1 = R*p_1;

  % Rotate p_2 to be in the xy-plane with p_1.
  p_2 = R*p_2;
  p2x = p_2(1); p2y = p_2(2); p2z = p_2(3);
  phi = atan2(p2z, p2x);
  p_2 = Ry(phi)*p_2; % TODO: WHY NO MINUS?
  R = Ry(phi)*R;

  p_3 = R*p_3;
  p_4 = -(p_1 + p_2 + p_3);

  % Put everything back into a row vector.
  p_1 = p_1'; p_2 = p_2'; p_3 = p_3'; p_4 = p_4';

  % Set the z components (and also y in case of carbon) to zero so
  % they all have exactly the same value rather than 0.0000 and
  % -0.0000, etc.
  p_1 = [p_1(1) 0 0];
  p_2 = [p_2(1:2) 0];

  out = [p_1 p_2 p_3 p_4];
end

function out = Rx(t)
  out = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t);]
end

function out = Ry(t)
  out = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t);]
end

function out = Rz(t)
  out = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1;]
end

% IonVelocities: Calculate velocities of ion given initial parameters
% Usage: IonVelocities(InitialConditions) where InitialConditions is a vector
% of the form [x1, y1, z1, x2, ..., z3].
function out = ionVelocities (initialConditions, masses, charges)
  options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
  [t,y] = ode45('hamiltonianDerivative', [0 1e-11], [initialConditions masses charges], options);
  out = y(size(t,1), 10:18);
end
