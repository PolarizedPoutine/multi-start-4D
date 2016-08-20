% Natural Extension Reference Frame algorithm
% Parsons et al., "Practical conversion from torsion space to Cartesian space
% for in silico protein synthesis", Journal of computational chemistry 26(10),
% 1063-1068 (2005).

% Z_i = [rC r_i tC t_i pC p_i]

function C = ZMatrixToCartesian(Z)
atoms = size(Z,1);
C = zeros(atoms,3);

if atoms >= 2
  r2 = Z(2,2);
  C(2,1) = r2;

  if atoms >= 3
    rC = uint64(Z(3,1));
    r3 = Z(3,2);
    theta = Z(3,4);

    % TODO: Assert what if rC == thetaC ?

    if rC == 1
      x2 = r3 * cosd(180 - theta);
      y2 = r3 * sind(180 - theta);
      C(3,1:2) = C(1,1:2) + [x2, y2];
    else % Must be 2 but TODO: Assert what if it's not 1 or 2!
      x2 = r3 * cosd(theta);
      y2 = r3 * sind(theta);
      C(3,1:2) = [x2, y2];
    end

    if atoms >= 4
      for i = 4:atoms
        rC = uint64(Z(i,1));
        r = Z(i,2);
        thetaC = uint64(Z(i,3));
        theta = Z(i,4);
        phiC = uint64(Z(i,5));
        phi = Z(i,6);

        x = r * cosd(theta);
        y = r * cosd(phi) * sind(theta);
        z = r * sind(phi) * sind(theta);

        r_rC = C(rC,:);
        r_thetaC = C(thetaC,:);
        r_phiC = C(phiC,:);

        ab = r_thetaC - r_phiC;
        bc = (r_rC - r_thetaC) / norm(r_rC - r_thetaC);
        n = cross(ab,bc); n = n / norm(n);
        ncbc = cross(n,bc);
        R = [bc' ncbc' n'];

        C(i,:) = (r_rC' + R*[-x; y; z;])';
      end
    end
end
end
