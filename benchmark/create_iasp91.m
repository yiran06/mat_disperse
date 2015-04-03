%
function create_iasp91

[laym,vpm,vsm,hkm,depm] = Read_iasp91;

write_model_rbh('iasp91',hkm,vpm,vsm,zeros(size(vsm))+3.0);

end





% Read in IASP91 model
% each layer, the velocity is uniform
function [laym,vpm,vsm,hkm,depm] = Read_iasp91

Depmax = 50; % km
dm = 2;        % gradient is discretized into dm thick layers


data = load('/home/yma/Datalib/EarthModel/iasp91.mdl');
depth0 = data(:,1);
vp0 = data(:,3);
vs0 = data(:,4);
n0 = length(depth0);

% Combine the layers with same velocity
depth = zeros(n0,1);
vp = zeros(n0,1);
vs = zeros(n0,1);

depth(1) = depth0(1);
vp(1) = vp0(1);
vs(1) = vs0(1);
n = 1;
for i=2:n0
    if vp0(i) == vp0(i-1) && vs0(i) == vs0(i-1)...
            && vp0(i)==vp0(i+1) && vs0(i) == vs0(i+1)
        continue;
    end
    n = n + 1;
    depth(n) = depth0(i);
    vp(n) = vp0(i);
    vs(n) = vs0(i);
    
    if depth0(i) >= Depmax
        break
    end
end

% Transform gradient model to homogeneous layer model
% use 1 km
vpm = zeros(n,1);
vsm = zeros(n,1);
hkm = zeros(n,1);

laym= 0;
% To the maximum of 1000 km
% loop through n depths
for i=1:n-1
    if vp(i) == vp(i+1) && vs(i) == vs(i+1)  % v is constant over the depth
        laym = laym + 1;
        vpm(laym) = vp(i);
        vsm(laym) = vs(i);
        hkm(laym) = depth(i+1) - depth(i);
    else
        if depth(i) == depth(i+1) % discontinuity
            continue
        end
        
        % Transfer the gradient to layered
        hk   = depth(i+1) - depth(i);
        devp = (vp(i+1) - vp(i))/hk;
        devs = (vs(i+1) - vs(i))/hk;
        
        nl   = ceil(hk/dm);     
        hkl  = [repmat(dm,nl-1,1);hk-(nl-1)*dm];
        depl = cumsum(hkl);    
        midl = depl - 0.5 * hkl;
        
        vpl  = vp(i) + devp * midl;
        vsl  = vs(i) + devs * midl;
        
        vpm(laym+1:laym+nl) = vpl;
        vsm(laym+1:laym+nl) = vsl;
        hkm(laym+1:laym+nl) = hkl;
        laym= laym+nl;
       
    end
end

vpm = vpm(1:laym);
vsm = vsm(1:laym);
hkm = hkm(1:laym);

tmp = cumsum(hkm);
depm= [0; tmp(1:end-1)];

end

