function [lamd,lamu] = updown_love(thk,cvs,om,k,z,layer)

% This function calculates the up-going and down-going matrices for the SH case
% Note that the function psv also calculates up-going and down-going matrices
% (called du) which are optimized for use in calculating the modified R/T
% coefficients. The matrices calculated in this function are used in calculating
% the displacement-stress vectors.

cvs2 = cvs(layer).^2;
depth = [0 ; cumsum(thk)];

k2 = k^2; om2 = om^2;
   
ks2 = om2/cvs2;
nus = sqrt(k2-ks2);
if imag(-i*nus) > 0;
    nus = -nus;
end   
   
 
lamd = exp(-nus*(z-depth(layer)));
lamu = 0;  
if layer <= length(thk)
    lamu = exp(-nus*(depth(layer+1)-z));
end