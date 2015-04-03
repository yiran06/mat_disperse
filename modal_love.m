function cr = modal_love(freq,thk,dns,cvs,crmin,crmax)

% This function calculates the modal love wave phase velocities in an elastic,
% vertically heterogeneous medium using search techniques.

% Establish global parameters
global TOL MAXROOT NUMINC

% Initialize a matrix to store modal phase velocities
cr = zeros(length(freq),MAXROOT);

% Loop through the frequencies
for j = 1:length(freq)
    numroot = 0;
    om = 2*pi*freq(j);
    
    % Establish the search parameters
    kmax = om/crmin;
    kmin = om/crmax;
    dk = (kmax - kmin)/NUMINC;
    
    % Establish the first and second points
    k1 = kmax;
    f1 = secular_love(k1,om,thk,dns,cvs);
    k2 = kmax - dk;
    f2 = secular_love(k2,om,thk,dns,cvs);
    
    % Establish an arbitrary high value for kold
    kold = 1.1*kmax;
    
    % Loop through the remaining points
    for m = 2:NUMINC-1
        k3 = kmax - m*dk;
        f3 = secular_love(k3,om,thk,dns,cvs);
        
        % Determine if a minimum is bracketed
        if (f2 < f1) & (f2 < f3)
            
            % Use golden search/parabolic interpolation to refine minimun
            [ktrial,ftrial] = fminbnd('secular_love',k3,k1,optimset('TolX',1e-12,'Display','off'),om,thk,dns,cvs);
            
            
            % Check to see if ktrial is a zero and different from the previous zero
            if (ftrial < TOL & abs((ktrial-kold)/kold) > 1e-2)
                numroot = numroot + 1;
                cr(j,numroot) = om/ktrial;
                kold = ktrial;
            end
        end
        
        % Break out of loop of maxroots is reached
        if numroot == MAXROOT
            break;
        end
        
        k1 = k2; f1 = f2;
        k2 = k3; f2 = f3;
        
    end
end