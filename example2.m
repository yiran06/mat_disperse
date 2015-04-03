% Input file
function example2
mdl_file = 'benchmark/CVM_1d.mdl';
fid = fopen(mdl_file);
for i=1:12
    fgetl(fid);
end
m = textscan(fid,'%f %f %f %f %d %d %d %d %d %d');
fclose(fid);

thk = m{1}; thk(end)=[];
dns = m{4};
vp  = m{2};
vs  = m{3};

% thk = [18.0 6.0 6.0];
% vp = [6 6.30 6.70 8.20];
% vs = [3.5 3.65 3.9 4.7];
% dns = [2.80 2.90 3.10 3.3];
% 
freq = 1./[1:20];


% Call mat_disperse.m to solve the eigenvalue problem and calculate phase
% velocities, displacement-stress functions, and surface wave displacements

% [vr,z,r,dvrvs,dvrvp,dvrrho] = mat_disperse(thk,dns,vp,vs,freq,'R');
% 
% subplot(211);
% for i = 1:length(freq)
%     plot(squeeze(abs(dvrvs(i,1,:))));
%     hold off;
% end

WAVE = 'L';
[vr,z,r,dvrvs,dvrrho] = mat_disperse(thk,dns,vp,vs,freq,WAVE);

if strcmp(WAVE,'R')
    test_with_rbh('RAYLEIGH',1./freq,vr,dvrvs,dvrrho);
elseif strcmp(WAVE,'L')
    test_with_rbh('LOVE',1./freq,vr,dvrvs,dvrrho);
end
 
end

function test_with_rbh(flag,T,vr,dvrvs,dvrrho)

if nargin < 3
    flag = 'LOVE';
    T = 1:2:20;
end


% read the test dispersion

% read dispersion
rfile = 'SREGN.ASC';
lfile = 'SLEGN.ASC';

[r_modes,r_period,r_C,r_U,r_ENERGY,~] = read_data(rfile);
[l_modes,l_period,l_C,l_U,l_ENERGY,~] = read_data(lfile);

r_id0 = find(r_modes == 0);
r_id1 = find(r_modes == 1);
l_id0= find(l_modes == 0);
l_id1= find(l_modes == 1);

figure(1);
plot(r_period(r_id0),r_C(r_id0),'ro-')
hold on;
plot(l_period(l_id0),l_C(l_id0),'bo-');
plot(T,vr,'co-');
title('Dispersion curve');
hold off;

figure(2);
for i=1:length(T)
    file = [num2str(T(i)),'.0.',flag,'.der.txt'];
    [dvp,dvs,dr,UT] = read_der(file);
    
    subplot(211)
    % plot dvs
    plot(dvs,'ko-');
    hold on;
    plot(squeeze(real(dvrvs(i,1,:))),'ro-');
    title(file)
    hold off;
     
    subplot(212)
    plot(dr,'ko-');
    hold on;
    plot(squeeze(real(dvrrho(i,1,:))),'ro-');
    hold off;
    
end
end

function [dvp,dvs,dr,UT] = read_der(file)

fid = fopen(['benchmark/',file]);
for i=1:3
    fgetl(fid);
end
line=fgetl(fid);
C = strsplit(line);
C(1) = [];
m = textscan(fid,repmat('%f ',1,length(C)));
fclose(fid);


[~,id] = ismember('DC/DB',C);
dvs = m{id};

[~,id] = ismember('DC/DR',C);
dr = m{id};

[flag,id] = ismember('DC/DR',C);
if flag
    dvp = m{id};
else
    dvp = NaN;
end

[flag,id] = ismember('UT',C);
if flag
    UT = m{id};
else
    UT = NaN;
end
end

function [modes,period,C,U,ENERGY,GAMMA] = read_data(file)

fid = fopen(['benchmark/',file]);
fgetl(fid);
if file(2)=='L'
m = textscan(fid,'%d %d %f %f %f %f %f %f');
else
m = textscan(fid,'%d %d %f %f %f %f %f %f %f');    
end
fclose(fid);

modes  = m{1};
period = m{3};
C = m{5};
U = m{6};
ENERGY = m{7};
GAMMA = m{8};
end