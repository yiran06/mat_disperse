function create_period
% 
% f = 0.1:0.1:1;
% T = 1./f; T = fliplr(T);
T = 1:1:20;
fid = fopen('period.txt','w+');
fprintf(fid,'%10.3f\n',T);
fclose(fid);