function multiStart4Ion(momenta, masses, charges, fOutFilenamePrefix, startingIndex, runs)
  nMomenta = size(momenta, 1);

  parfor i = startingIndex:nMomenta
    % fOutFilename = strcat(fOutFilenamePrefix, '_G', sprintf('%05d',i), '.log');
    % fOut = fopen(fOutFilename, 'a');

    p = momenta(i,:)
    p = removeCOMMotion4Ion(p, masses)
    p = rotateMomentum4Ion(p)

    pGoal = p;
    residualNormObjective = @(g)residualNorm(g, pGoal, masses, charges);

    lowerBounds = [0 0 0 0 0 0];
    upperBounds = [1000 1000 1000 180 180 180];

    initialGeometry = [106 120 106 172 175 10];

    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', ...
      'iter', 'MaxFunEvals', 3000);
    problem = createOptimProblem('fmincon', 'objective', residualNormObjective, ...
      'lb', lowerBounds, 'ub', upperBounds, 'x0', initialGeometry, 'options', ...
      options);
    ms = MultiStart('UseParallel', 'always', 'Display', 'iter', 'StartPointsToRun', 'bounds');
    [g, fval, exitflag, output, solutions] = run(ms, problem, runs);

    fprintf('G%05d DONE @ %s.\n', i, datestr(now));
    %fprintf('G%05d Best geometry found: (%.2f pm, %.2f pm, %.2f deg) with log residual norm %.2f and exit flag %d.\n', i, g(1), g(2), g(3), fval, exitflag);
    fprintf('G%05d Solver: funcCountNumber:       %d\n', i, output.funcCount);
    fprintf('G%05d         localSolverIncomplete: %d\n', i, output.localSolverIncomplete);
    fprintf('G%05d         localSolverNoSolution: %d\n', i, output.localSolverNoSolution);
    fprintf('G%05d         localSolverSuccess:    %d\n', i, output.localSolverSuccess);
    fprintf('G%05d         localSolverTotal:      %d\n', i, output.localSolverTotal);

    % We don't always get back 'runs' solutions so we count how many we found.
    numSolutionsFound = size([solutions.Fval], 2);

    % We put each distinct solution we found into a row vector and  print them all.
    mostLikelyGeometries = [i*ones(numSolutionsFound, 1) reshape([solutions.X], 6, numSolutionsFound)' [solutions.Fval]' [solutions.Exitflag]'];

    %fprintf('G%05d Writing %s.\n', i, fOutFilename);
    for j = 1:numSolutionsFound
      fprintf('G%05d G %d\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%2.2f\t%d\n', i, mostLikelyGeometries(j,:));
      %fprintf(fOut, '%d\t%3.6f\t%3.6f\t%3.6f\t%2.2f\t%d\n', mostLikelyGeometries(j,:));
    end

    % fprintf('\n');
    % fprintf(fOut,'\n');
    % fclose(fOut);
  end
end

% This is our objective or fitness function for the multi start algorithm. It
% takes a vector g with our molecule's configuration (r_12, r_23, theta) and
% simulates a Coulomb explosion for it. It then compares the asymptotic
% momentum produced with the goal momentum we are attempting to recreate and
% returns the log10 of the norm of the difference between the two momentum
% squared.
function rn = residualNorm(g, pGoal, masses, charges)
  g = [1e-12*g(1:3) g(4:6)];
  p = simulateMomentum4Ion(g, masses, charges, false);
  rn = log10(norm(pGoal - p)^2);
end
