function [ur,uy] = green(freq,vr,group,energy,z,r,offsets,Fx,Fy,Fz,Phi,s_depth,r_depth)

% This function calculates the Green's function of the horizontal and vertical
% Rayleigh wave displacement field using modal superposition.

% Copyright 1999 by Glenn J. Rix and Carlo G. Lai

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Initialize matrices for the horizontal and vertical displacements
ur = zeros(length(freq),length(offsets));
uy = zeros(length(freq),length(offsets));

% Define the vector of circular frequencies
om = 2*pi*freq;

% Loop through the frequencies
for j = 1:length(freq)
   
   % Calculate the position of the source and receiver
   index = find(z(:,j) >= s_depth);
   s_index = index(1);
   index = find(z(:,j) >= r_depth);
   r_index = index(1);
   
   % Determine the number of modes at each frequency
   index = find(vr(j,:));
   
   % Calculate the wavenumbers
   k = om(j)./shiftdim(vr(j,index),1);

% Assign displacement vectors, velocities, and energy integrals to local variables
   r1s = shiftdim(r(j,index,s_index,1),1);
   r1r = shiftdim(r(j,index,r_index,1),1);
   r2s = shiftdim(r(j,index,s_index,2),1);
   r2r = shiftdim(r(j,index,r_index,2),1);
   V = shiftdim(vr(j,index),1);
   U = shiftdim(group(j,index),1);
   I1 = shiftdim(energy(j,index),1);
   
   % Loop through the offsets
   for n = 1:length(offsets)
      a = (Fz*r2s + i*(Fx*cos(Phi)+Fy*sin(Phi))*r1s)./(8*V.*U.*I1)./...
         sqrt((pi*k*offsets(n))/2).*exp(-i*k*offsets(n)).*r1r*exp(-i*pi/4);
      b = (Fz*r2s + i*(Fx*cos(Phi)+Fy*sin(Phi))*r1s)./(8*V.*U.*I1)./...
         sqrt((pi*k*offsets(n))/2).*exp(-i*k*offsets(n)).*r2r*exp(i*pi/4);
      ur(j,n) = sum(a);
      uy(j,n) = sum(b);
   end
end   
