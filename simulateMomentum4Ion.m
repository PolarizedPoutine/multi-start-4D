function out = simulateMomentum4Ion(g, masses, charges)
  p_0 = zeros(1,12);
  p = timeEvolveMomenta4Ion([g p_0], masses, charges);
  out = rotateMomentum4Ion(p);
end

function out = timeEvolveMomenta4Ion(initialConditions, masses, charges)
  options = odeset('AbsTol', 1e-27, 'RelTol', 1e-6, 'InitialStep', 1e-18);
  [t,y] = ode45('hamiltonsEquations4Ion', [0 1e-11], [initialConditions masses charges], options);
  out = y(size(t,1), 13:24);
end
