function create_model
z = 0:1:40;
h1= 1;
name = 'CVM_1d.mdl';

vp_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vp_basin.grd';
vs_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vs_basin.grd';
rho_file= '/export/nobackup/yma/LASSIE/CVMh_ref/rho_basin.grd';


% vp_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vp_1d.grd';
% vs_file = '/export/nobackup/yma/LASSIE/CVMh_ref/vs_1d.grd';
% rho_file= '/export/nobackup/yma/LASSIE/CVMh_ref/rho_1d.grd';




vp_data = load(vp_file);
vs_data = load(vs_file);
rho_data= load(rho_file);

vp = interpp(vp_data(:,1),vp_data(:,2),z);
vs = interpp(vs_data(:,1),vs_data(:,2),z);
rho= interpp(rho_data(:,1),rho_data(:,2),z);

[h,vp2,vs2,rho2] = convert_z_h(z,vp,vs * 80/100,rho,h1);

write_model_rbh(name,h,vp2,vs2,rho2);

right = cumsum(h); right_vp = vp2; right_vs = vs2;
left  = right - h(1); left_vp = vp2; left_vs = vs2;
all = zeros(1,length(left) * 2);
all_vp = all;
all_vs = all;
all(1:2:end) = left;
all(2:2:end) = right;
all_vp(1:2:end) = left_vp;
all_vp(2:2:end) = right_vp;
all_vs(1:2:end) = left_vs;
all_vs(2:2:end) = right_vs;

plot(all_vp,all,'r-');
hold on;
plot(all_vs,all,'b-');
hold off;
legend('vp','vs','location','southwest');
xlabel('km/s');ylabel('Depth (km)');
grid on;
set(gca,'YDir','reverse');



% convert depth-v to h-v
function [h,vp2,vs2,rho2] = convert_z_h(z,vp,vs,rho,h1)

dep1 = z(1):h1:z(end)-h1;
dep2 = dep1+h1;
vp1  = interpp(z,vp,dep1);  vp2  = interpp(z,vp,dep2);
vs1  = interpp(z,vs,dep1);  vs2  = interpp(z,vs,dep2);
rho1 = interpp(z,rho,dep1); rho2 = interpp(z,rho,dep2);

h = dep2 - dep1;
vp2  = (vp1 + vp2)/2;
vs2  = (vs1 + vs2)/2;
rho2 = (rho1 + rho2)/2;


% write the model
function write_model_rbh(name,h,vp,vs,rho)
fid = fopen(name,'w+');
fprintf(fid,'MODEL\nTEST MODEL\nISOTROPIC\nKGS\nFLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\nH(KM)  VP  VS  RHO  QP   QS  ETAP  ETAS  FREFP  FREFS \n');

for i=1:length(h)
   fprintf(fid,'%10.3f %10.3f %10.3f %10.3f 1400 600 0 0 1 1\n',h(i),vp(i),vs(i),rho(i));
end
fprintf(fid,'%10.3f %10.3f %10.3f %10.3f 1400 600 0 0 1 1\n',0,vp(end),vs(end),rho(end));

fclose(fid);



function y1 = interpp(x0,y0,x1)
y1 = interp1(x0,y0,x1);

