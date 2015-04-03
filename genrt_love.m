function [Td,Rd] = genrt_love(td, tu, rd, ru)

N = length(td);

Td = zeros(N,1);
Rd = zeros(N,1);

% Calculate the Td and Rd matrices for the Nth layer
Td(N) = td(N);
Rd(N) = rd(N);

% Loop through the first N-1 layers in reverse order
for j = N-1:-1:1
   Td(j) = td(j)/(1 - ru(j)*Rd(j+1));
   Rd(j) = rd(j) + tu(j)*Rd(j+1)*Td(j);
end