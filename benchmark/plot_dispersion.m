function plot_dispersion

% read dispersion
rfile = 'SREGN.ASC';
lfile = 'SLEGN.ASC';

[r_modes,r_period,r_C,r_U,r_ENERGY,~] = read_data(rfile);
[l_modes,l_period,l_C,l_U,l_ENERGY,~] = read_data(lfile);

r_id0 = find(r_modes == 0);
r_id1 = find(r_modes == 1);
l_id0= find(l_modes == 0);
l_id1= find(l_modes == 1);


% rayleigh in red, love in blue
% 0 mode in -, 1 mode in .-
% phase velocity
figure
set(gcf,'position', [21 32 500 1065],'PaperPositionMode','auto');
subplot(211)
plot(r_period(r_id0),r_C(r_id0),'ro-')
hold on;
plot(l_period(l_id0),l_C(l_id0),'bo-');
plot(r_period(r_id1),r_C(r_id1),'r*-')
plot(l_period(l_id1),l_C(l_id1),'b*-');
xlabel('Period (s)');
ylabel('(km/s)');
xlim([0 10]);
% ylim([0 2]);
grid on;box on;
legend('R0','L0','R1','L1','Location','northwest');
title('Phase Velocity');
hold off;

subplot(212)
plot(r_period(r_id0),r_U(r_id0),'ro-')
hold on;
plot(l_period(l_id0),l_U(l_id0),'bo-');
plot(r_period(r_id1),r_U(r_id1),'r*-')
plot(l_period(l_id1),l_U(l_id1),'b*-');
xlabel('Period (s)');
ylabel('(km/s)');
xlim([0 10]);
% ylim([0 2]);
grid on;box on;
legend('R0','L0','R1','L1','Location','northwest');
title('Group Velocity');
hold off;

end

function [modes,period,C,U,ENERGY,GAMMA] = read_data(file)

fid = fopen(file);
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