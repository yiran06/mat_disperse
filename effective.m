function [vre,dvrevs,dvrevp] = effective(freq,vr,group,energy,z,r,dvrvs,dvrvp,offsets,s_depth,r_depth)

% This function calculates the effective vertical phase velocities resulting
% from the superposition of the complex-valued modal phase velocities. 
% Only the case of a vertical point load (i.e., Fxx = Fyy = 0) is considered. The
% algorithm is based on Eqs. 3.46 and 3.57 of Lai (1998).

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

% Determine the number of layers
[m n p] = size(dvrvs);

% Initialize matrices for the effective vertical velocities and their derivatives
vre = zeros(length(freq),length(offsets));
dvrevs = zeros(length(freq),length(offsets),p);
dvrevp = zeros(length(freq),length(offsets),p);

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
   m = length(index);
   
   % Initialize matrices
   A = zeros(m); B = zeros(m); Lambda = zeros(m); Sigma = zeros(m); Chi = zeros(m); Eta = zeros(m);
   E1 = zeros(m,m,p); E2 = zeros(m,m,p);
   
   % Assign displacement vectors, velocities, and energy integrals to local variables
   r2r = squeeze(r(j,index,r_index,2));
   r2s = squeeze(r(j,index,s_index,2));
   V = squeeze(vr(j,index));
   U = squeeze(group(j,index));
   I1 = squeeze(energy(j,index));
   
   % Loop through the offsets
   for n = 1:length(offsets)
      
      % Terms in Eqs. 3.46
      for m1 = 1:length(index)
         for m2 = 1:length(index)
            Psiy = (r2r(m1)*r2r(m2)*r2s(m1)*r2s(m2))/(U(m1)*I1(m1)*U(m2)*I1(m2));
            B(m1,m2) = Psiy*cos(om(j)*offsets(n)*(1/V(m1)-1/V(m2)))/sqrt(V(m1)*V(m2));
            A(m1,m2) = B(m1,m2)*(1/V(m1)+1/V(m2));
            C = sin(om(j)*offsets(n)*(1/V(m1)-1/V(m2)));
            D = cos(om(j)*offsets(n)*(1/V(m1)-1/V(m2)));
            Lambda(m1,m2) = Psiy*(2*om(j)*offsets(n)*C - V(m1)*D)*...
               V(m2)^2/(2*(V(m1)*V(m2))^2*sqrt(V(m1)*V(m2)));
            Sigma(m1,m2) = Psiy*(2*om(j)*offsets(n)*C + V(m2)*D)*...
               V(m1)^2/(2*(V(m1)*V(m2))^2*sqrt(V(m1)*V(m2)));
            Chi(m1,m2) = Psiy*(2*(om(j)*offsets(n)*V(m2)*C - V(m1)*V(m2)*...
               D + om(j)*offsets(n)*V(m1)*C) - V(m1)*D*(V(m1)+V(m2)))*...
               V(m2)^2/(2*(V(m1)*V(m2))^3*sqrt(V(m1)*V(m2)));
            Eta(m1,m2) = Psiy*(2*(om(j)*offsets(n)*V(m2)*C + V(m1)*V(m2)*...
               D + om(j)*offsets(n)*V(m1)*C) + V(m2)*D*(V(m1)+V(m2)))*...
               V(m1)^2/(2*(V(m1)*V(m2))^3*sqrt(V(m1)*V(m2)));
         end
      end
      
      % Equation B.12
      A1 = sum(sum(A));
      B1 = sum(sum(B));
      
      % Equation 3.46
      vre(j,n) = 2*B1/A1;
      
      % Equation B.16
      FF = 2*(A1*Lambda - B1*Chi)/A1^2;
      GG = 2*(A1*Sigma - B1*Eta)/A1^2;
      
      for m1 = 1:length(index)
         for m2 = 1:length(index)
            E1(m1,m2,:) = FF(m1,m2)*dvrvs(j,m1,:) - GG(m1,m2)*dvrvs(j,m2,:);
            E2(m1,m2,:) = FF(m1,m2)*dvrvp(j,m1,:) - GG(m1,m2)*dvrvp(j,m2,:);
         end
      end
      
      dvrevs(j,n,:) = shiftdim(sum(sum(E1,1),2),2);
      dvrevp(j,n,:) = shiftdim(sum(sum(E2,1),2),2);
      
   end
end   
