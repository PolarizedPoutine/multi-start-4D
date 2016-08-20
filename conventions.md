# Conventions
A list of conventions I follow throughout the project.

## Common variables.
Unless otherwise stated, positions are labeled r and momentum p.

## Units
All functions should accept input and produce output in SI units. Only
exceptions are masses (in amu) and atomic charges (in units of the elementary
charge e) for clarity. Functions, especially mutliStart4Ion, may convert bonds
lengths to picometers so that all geometrical parameters are on the same order
of magnitude.

## Momentum vectors
The momentum vectors will always be in the order CCHH.

### Momentum vector orientation
First atom (C) along the x-axis. Second atom (C) in the xy-plane so both C's
form a plane. Third atom (H) is rotated using the same rotation that brought the
first two atoms into the xy-plane so that it sticks out of the xy-plane. Fourth
atom (H) is calculated using conservation of momentum.

## Geometries
May be stored as a Z-matrix or as a 6-element row vector.
TODO: - [] When do you use Z-matrix? When to use the 6-element row vector?

### Positions of the atoms
Convert to Z-matrix then to Cartesian? Because this way we get a unique
Cartesian representation for each geometry described.

## Z-matrix
For 4-ion acetylene the [Z-matrix] (https://en.wikipedia.org/wiki/Z-matrix_(chemistry))
has 6 parameters. I will always call these individual parameters r12, r23, r34,
theta13, theta24, and phi.

The Z-matrix for acetylene at equilibrium requires
the use of dummy atoms due to the 180° bond angles but we are looking at the
dynamics of acetylene so we cannot assume a 180° bond angle so we don't need
to use dummy atoms. Ummm wait a second, my Z-matrix works with 180° bond angles.
I should check this out.
