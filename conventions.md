# Conventions
A list of conventions I follow in this code.

## Units
All functions should accept input and produce output in SI units. Only
exceptions are masses (in amu) and atomic charges (in units of the elementary
charge e) for clarity. Functions, especially mutliStart4Ion, may convert bonds
lengths to picometers so that all geometrical parameters are on the same order
of magnitude.

## Momentum vectors
The momentum vectors will always be in the order CCHH.

## Geometries
May be stored as a Z-matrix or as a 6-element row vector.

## Z-matrix
For 4-ion acetylene the [Z-matrix] (https://en.wikipedia.org/wiki/Z-matrix_(chemistry))
has 6 parameters. I will always call these individual parameters r12, r23, r34,
theta13, theta24, and phi.
