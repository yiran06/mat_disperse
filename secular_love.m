function d = secular_love(k,om,thk,dns,cvs)

% epsilon = 0.0001;
% while any(abs(om/k-cvs)<epsilon) | any(abs(om/k-cvp)<epsilon)
%    k = k * (1+epsilon);
% end   

[e11,e12,e21,e22,du] = sh(thk,dns,cvs,om,k);
[td,tu,rd,ru] = modrt_love(e11,e12,e21,e22,du);
[Td, Rd] = genrt_love(td, tu, rd, ru);

Ru_0 = -e22(1)*du(1)/e21(1);
d = abs( 1-Ru_0*Rd(1));
