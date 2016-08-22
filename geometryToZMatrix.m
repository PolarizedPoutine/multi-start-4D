function out = geometryToZMatrix(g)
  r_H1C1 = g(1); r_C1C2 = g(2); r_C2H2 = g(3);
  theta_H1C2 = g(4); theta_C1H2 = g(5); phi_H1H2 = g(6);

  out = [0 0      0 0          0 0;        ... % C1
         1 r_C1C2 0 0          0 0;        ... % C2
         1 r_H1C1 2 theta_H1C2 0 0;        ... % H1
         2 r_C2H2 1 theta_C1H2 3 phi_H1H2; ];  % H2
end
