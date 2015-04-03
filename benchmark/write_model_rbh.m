% write the model
function write_model_rbh(name,h,vp,vs,rho)
fid = fopen(name,'w+');
fprintf(fid,'MODEL\nTEST MODEL\nISOTROPIC\nKGS\nFLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\nH(KM)  VP  VS  RHO  QP   QS  ETAP  ETAS  FREFP  FREFS \n');

for i=1:length(h)
   fprintf(fid,'%10.3f %10.3f %10.3f %10.3f 1400 600 0 0 1 1\n',h(i),vp(i),vs(i),rho(i));
end
fprintf(fid,'%10.3f %10.3f %10.3f %10.3f 1400 600 0 0 1 1\n',0,vp(end),vs(end),rho(end));

fclose(fid);
