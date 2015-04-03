function [I1,I2,I3,U,zdvrvs,zdvrrho,dvrvs,dvrrho] = partial_love(freq,vr,z,r,dr,thk,dns,cvs)

global NUMPOINTS MAXROOT

% Initialize vectors for material properties
vs = zeros(NUMPOINTS,1);
shear = zeros(NUMPOINTS,1);
rho = zeros(NUMPOINTS,1);

% Initialize matrices for the energy integrals, group velocity, and partial derivatives
I1 = zeros(length(freq),MAXROOT);
I2 = zeros(length(freq),MAXROOT);
I3 = zeros(length(freq),MAXROOT);
U = zeros(length(freq),MAXROOT);
zdvrvs = zeros(length(freq),MAXROOT,NUMPOINTS); zdvrrho = zdvrvs;
dvrvs = zeros(length(freq),MAXROOT,length(dns)); dvrrho = dvrvs;

% Define the vector of circular frequencies
om = 2*pi*freq;

% Loop through the frequencies
for j = 1:length(freq)
   
   % Loop through the vector of depths to assign vectors of material properties
   for n = 1:NUMPOINTS
      
      % Determine the layer corresponding to the current depth
      index1 = find(z(n,j) <= [cumsum(thk) ; z(NUMPOINTS,j)]);
      layer = index1(1);
      
      % Assign layer properties to vectors
      vs(n) = cvs(layer);
      shear(n) = dns(layer)*cvs(layer)*cvs(layer);
      rho(n) = dns(layer);
   end      
   
   % Loop through the modes at each frequency
   index2 = find(vr(j,:));
   for m = 1:length(index2)
      
      % Calculate the wavenumber
      k = om(j)/vr(j,index2(m));
      
      % Assign the displacement vectors and their derivatives to local variables
      r1 = squeeze(r(j,m,:,1));
      dr1 = squeeze(dr(j,m,:,1));
      
      % Calculate the first energy integral
      integrand = rho.*(r1.^2);
      I1(j,m) = 0.5*trapz(z(:,j),integrand);
      
      % Calculate the second energy integral
      integrand = shear.*(r1.^2);
      I2(j,m) = 0.5*trapz(z(:,j),integrand);
      
      % Calculate the third energy integral
      integrand = shear.*(dr1.^2);
      I3(j,m) = 0.5*trapz(z(:,j),integrand);
      
      % Calculate the group velocity
      U(j,m) = I2(j,m)/I1(j,m)/vr(j,index2(m));
      
      % Calculate the partial derivatives at each individual depth
      zdvrvs(j,m,:)  = vr(j,index2(m)).*rho.*vs.*( k^2*r1.^2 + dr1.^2 )/(2*k^2*I2(j,m));
      zdvrrho(j,m,:) = vr(j,index2(m)).*( (k^2*r1.^2 + dr1.^2).*vs.^2 -om(j)^2.*r1.^2 )/(4*k^2*I2(j,m));
      
      % add the layer boundary points
      znew    = [z(:,j); cumsum(thk)]; znew = unique(znew);
      dvsnew  = interp1(z(:,j), squeeze(zdvrvs(j,m,:)), znew);
      drhonew = interp1(z(:,j), squeeze(zdvrrho(j,m,:)), znew);
      
%       plot(z(:,j),squeeze(zdvrvs(j,m,:)),'ko-');
%       hold on;
%       plot(znew,dvsnew,'r-');
%       xlim([0 10]);
%       hold off;
      
      % Calculate the partial derivatives for each layer by integrating over the
      % thickness of the layer
      depth = [0 ; cumsum(thk) ; z(NUMPOINTS,j)];
      for n = 1:length(dns)
         index3 = find(znew >= depth(n) & znew <= depth(n+1));
         if length(index3) < 5
            disp(['Partial derivatives at ',num2str(freq(j)),' Hz for Layer ',num2str(n), ' may be incorrect.'])
         end   
         dvrvs(j,m,n)  = trapz(znew(index3),dvsnew(index3),1);
         dvrrho(j,m,n) = trapz(znew(index3),drhonew(index3),1);
      end
   end
end

end