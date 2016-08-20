function out = geometryToZMatrix(g)
  r12 = g(1); r23 = g(2); r34 = g(3);
  theta13 = g(4); theta24 = g(5); phi = g(6);

  % TODO: Pick Z-matrix!
  out = [0 0   0 0       0 0;   ... % H -
         n r12 0 0       0 0;   ... % - C -
         n r23 m theta13 0 0;   ... % - C -
         n r34 l theta24 k phi; ]   % - H
end
