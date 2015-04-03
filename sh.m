function [e11,e12,e21,e22,du] = sh(thk,dns,cvs,om,k)

cvs2 = cvs.^2;
mu = dns.*cvs2;

k2 = k.^2; om2 = om.^2;
ks2 = om2./cvs2;
nus = sqrt(k2-ks2);
index = find(imag(-i*nus) > 0);
nus(index) = -nus(index);

N = length(cvs);

e11= zeros(N,1)+1;
e12= zeros(N,1)+1;
e21= -mu.*nus;
e22= -e21;
du = exp(-nus(1:length(thk)).*thk);
