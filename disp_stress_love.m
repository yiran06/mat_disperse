function [z,r,dr] = disp_stress_love(freq,vr,thk,dns,cvs)

% Establish global parameters
global NUMPOINTS LAMBDA MAXROOT

% Calculate the vector of circular frequencies
om = 2*pi*freq;

% Determine the number of layers not including the half space
N = length(thk);

% Calulate the maximum depth for determining the displacement-stress vectors
lambda_max = max(sum(thk)+thk(end),LAMBDA*(real(vr(:,1)).^2 + imag(vr(:,1)).^2)./real(vr(:,1))./freq);
lambda_max = round(lambda_max);
% Initiate the depth and displacement-stress vectors and their numerical derivatives
z = zeros(NUMPOINTS,length(freq));
r = zeros(length(freq),MAXROOT,NUMPOINTS,2);
dr = zeros(length(freq),MAXROOT,NUMPOINTS,1);



% Loop through the frequencies
for j = 1:length(freq)
    
    % Create a vector of depths
    z(:,j) = linspace(0,lambda_max(j),NUMPOINTS)';
    step = floor(0.1/(z(2,j)-z(1,j)));
    n_try = [1:step:NUMPOINTS-1 NUMPOINTS];
    
    % Loop through the modes at each frequency
    index1 = find(vr(j,:));
    for m = 1:length(index1)
        
        % Calculate the wavenumber and load vector
        k = om(j)/vr(j,index1(m));
        
        % Check to see if the phase velocity is equal to the shear wave velocity
        % or compression wave velocity of one of the layers
        epsilon = 0.0001;
        while any(abs(om(j)/k-cvs)<epsilon) %| any(abs(om(j)/k-cvp)<epsilon)
            k = k * (1+epsilon);
        end
        
        % Calculate the SH element matrices for each layer and generalized R/T matrices
        [e11,e12,e21,e22,du] = sh(thk,dns,cvs,om(j),k);
        [td,tu,rd,ru] = modrt_love(e11,e12,e21,e22,du);
        [Td, Rd] = genrt_love(td, tu, rd, ru);
        
        % Initialize the Cd and Cu matrices
        cd = zeros(N+1,1);
        cu = zeros(N+1,1);
       
        cd(1) = 1;
        
        % Calculate Cd and Cu for the remaining layers
        for n = 1:N
            cu(n) = Rd(n)*cd(n);
            cd(n+1) = Td(n)*cd(n);
        end
        
        % Loop through the vector of depths
        for n = n_try
            % Determine the layer corresponding to the current depth
            index2 = find(z(n,j) <= [cumsum(thk) ; z(NUMPOINTS,j)]);
            layer = index2(1);
            
            % Calculate the up-going and down-going matrices for this depth
            [lamd,lamu] = updown_love(thk,cvs,om(j),k,z(n,j),layer);
            
            % Calculate the displacement-stress vector
            r(j,m,n,:) = [e11(layer) e12(layer) ; e21(layer) e22(layer)] * ...
                [lamd 0 ; 0 lamu] * [cd(layer) ; cu(layer)];
        end
        
        r1_try = abs(squeeze(r(j,m,n_try,1)));
        id = max(find(abs(diff(r1_try)) > 1e-6));
        
        for n = 1:n_try(id)
            % Determine the layer corresponding to the current depth
            index2 = find(z(n,j) <= [cumsum(thk) ; z(NUMPOINTS,j)]);
            layer = index2(1);
            
            % Calculate the up-going and down-going matrices for this depth
            [lamd,lamu] = updown_love(thk,cvs,om(j),k,z(n,j),layer);
            
            % Calculate the displacement-stress vector
            r(j,m,n,:) = [e11(layer) e12(layer) ; e21(layer) e22(layer)] * ...
                [lamd 0 ; 0 lamu] * [cd(layer) ; cu(layer)];
        end
        
        if n_try(id)~=NUMPOINTS
            r(j,m,n_try(id+1):NUMPOINTS,1) = interp1(  n_try(id+1:end),squeeze(r(j,m,n_try(id+1:end),1)),n_try(id+1):NUMPOINTS );
            r(j,m,n_try(id+1):NUMPOINTS,2) = interp1(  n_try(id+1:end),squeeze(r(j,m,n_try(id+1:end),2)),n_try(id+1):NUMPOINTS );
        end
        
        % Calculate the numerical derivative of the displacement-stress vectors
        % Note that only dr1 and dr2 are needed later. dr3 and dr4 are not calculated.
        dr(j,m,:,:) = gradient(squeeze(r(j,m,:,1)),z(:,j));
        
    end
end
