function plotTrajectory4Ion(t, y, masses, charges)

  amu = 1.66053886e-27;
  e = 1.60217646e-19;

  mC = amu*masses(1); mH = amu*masses(3);
  qC1 = charges(1); qC2 = charges(2); qH1 = charges(3); qH2 = charges(4);

  C1 = y(:,1:3);
  C2 = y(:,4:6);
  H1 = y(:,7:9);
  H2 = y(:,10:12);

  pC1 = y(:,13:15);
  pC2 = y(:,16:18);
  pH1 = y(:,19:21);
  pH2 = y(:,22:24);

  EC1 = (sqrt(sum(pC1.^2, 2)).^2 ./ (2*mC)) / e;
  EC2 = (sqrt(sum(pC2.^2, 2)).^2 ./ (2*mC)) / e;
  EH1 = (sqrt(sum(pH1.^2, 2)).^2 ./ (2*mH)) / e;
  EH2 = (sqrt(sum(pH2.^2, 2)).^2 ./ (2*mH)) / e;

  subplot(2,2,1);
  plot3(C1(:,1), C1(:,2), C1(:,3), '-o', C2(:,1), C2(:,2), C2(:,3), '-o', H1(:,1), H1(:,2), H1(:,3), '-o', H2(:,1), H2(:,2), H2(:,3), '-o');
  legend('C1', 'C2', 'H1', 'H2', 'Location', 'SouthEast');
  axis([-5e-10 5e-10 -5e-10 5e-10 -1e-10 1e-10]);
  title('Ion position (short timescale)');
  xlabel('X (m)');
  ylabel('Y (m)');
  grid on;

  subplot(2,2,2);
  plot3(C1(:,1), C1(:,2), C1(:,3), '-o', C2(:,1), C2(:,2), C2(:,3), '-o', H1(:,1), H1(:,2), H1(:,3), '-o', H2(:,1), H2(:,2), H2(:,3), '-o');
  legend('C1', 'C2', 'H1', 'H2', 'Location', 'SouthEast');
  title('Ion position (long timescale)');
  xlabel('X (m)');
  ylabel('Y (m)');
  grid on;

  subplot(2,2,3);
  plot3(pC1(:,1), pC1(:,2), pC1(:,3), '-o', pC2(:,1), pC2(:,2), pC2(:,3), '-o', pH1(:,1), pH1(:,2), pH1(:,3), '-o', pH2(:,1), pH2(:,2), pH2(:,3), '-o');
  legend('C1', 'C2', 'H1', 'H2', 'Location', 'SouthEast');
  title('Ion momentum');
  xlabel('X momentum (kg m/s)');
  ylabel('Y momentum (kg m/s)');
  grid on;

  subplot(2,2,4);
  plot(t, EC1, t, EC2, t, EH1, t, EH2, 'LineWidth', 1);
  legend('C1', 'C2', 'H1', 'H2', 'Location', 'SouthEast');
  title('Ion kinetic energy');
  xlabel('Time (s)');
  ylabel('Energy (eV)');
  xlim([0 1e-13]);
  grid on;

end
