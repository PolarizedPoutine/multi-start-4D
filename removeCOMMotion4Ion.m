% removeCOMMOtion takes a momentum triple, [p_1 p_2 p_3] and returns the triple
% in the same order with the center of mass motion removed.
% masses is the mass of each of the atoms (O, then C, then S) in amu.
function out = removeCOMMotion(momentum, masses)
  % We split our momentum triple into X,Y,Z components.
  p_x = momentum(1:3:10);
  p_y = momentum(2:3:11);
  p_z = momentum(3:3:12);

  masses = [15.9994 12.0107 32.065];  % OCS
  massSum = sum(masses);

  % We eliminate COM motion to put the particles in the COM frame of motion.
  p_x = p_x - sum(p_x) .* (masses)/massSum;
  p_y = p_y - sum(p_y) .* (masses)/massSum;
  p_z = p_z - sum(p_y) .* (masses)/massSum;

  % Putting the vectors back into the original format.
  momentum(1:12) = reshape([p_X; p_Y; p_Z], 1, 12);
  out = momentum;
end
