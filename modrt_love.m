function [td,tu,rd,ru] = modrt_love(e11,e12,e21,e22,du)

% This function calculates the modified R/T coefficients
N = length(du);

% Loop through the first N-1 layers
X = zeros(2,2,N);
for j = 1:N-1
    A = [e11(j+1) -e12(j);   e21(j+1) -e22(j)];
    B = [e11(j) -e12(j+1); e21(j) -e22(j+1)];
    L = [du(j) 0; 0 du(j+1)];
    X(:,:,j) = A\(B*L);
end

% Calculate the Nth layer
A = [e11(N+1) -e12(N); e21(N+1) -e22(N)];
B = [e11(N)*du(N); e21(N)*du(N)];
X(:,1,N) = A\B;

% Extract R/T submatrices
td = X(1,1,:);
ru = X(1,2,:);
rd = X(2,1,:);
tu = X(2,2,:);
