% Input file for mat_disperse.m

% Define vectors or layer thickness, mass density, shear wave velocity
% and compression wave velocity
thk = [5.0 10.0 10.0];
dns = [1.7 1.8 1.8 1.8];
vs = [200 300 400 500];
vp = [400 600 800 1000];

% Define a vector of frequencies (in Hz)
freq = linspace(1,10,10);

% Define a vector of offsets from the source
offsets = linspace(5,100,20);

% Call mat_disperse.m to solve the eigenvalue problem and calculate phase
% velocities, displacement-stress functions, and surface wave displacements
[vr,z,r,dvrvs,vre,dvrevs,ur,uy] = mat_disperse(thk,dns,vp,vs,freq,offsets);
