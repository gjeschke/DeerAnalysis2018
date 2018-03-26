function mk_ftikreg_input_std_alone(t,ff,regpar)
%
% Convert experimental data format (deer_fit input) to internal formats
%
% fname     file name, output of hermit interpolation

regmodstr='212 REGMOD'; % was 213 before May 17th 2005
lambdastr=[num2str(regpar) ' LAMBDA'];

% Range of analysis
%
rmin=1.5; 
rmax=8;
tdip=t;
dipevo=ff;

dospath=which('ftikreg_r_new.exe');
dospath=dospath(1:length(dospath)-length('ftikreg_r_new.exe'));

data=[tdip' dipevo'];
save([dospath 'ftikreg.dat'],'data','-ascii');
nexp=length(tdip);
ns=round(nexp/2);
ns=2*ns;
% Prepare FTIKREG_R parameter file
wfile=fopen([dospath 'ftikreg_r.par'],'w+');
fprintf(wfile,'%u%s\n',ns,' NS');
fprintf(wfile,'%3.1f%s\n',rmin,' SMIN');
fprintf(wfile,'%3.1f%s\n',rmax,' SMAX');
fprintf(wfile,'%s\n','1 DISMOD');
fprintf(wfile,'%s\n','0 M');
fprintf(wfile,'%u%s\n',nexp,' N');
fprintf(wfile,'%s\n','1 NE');
fprintf(wfile,'%s\n','22 ERRMOD');
fprintf(wfile,'%s\n','0.03 ERROR');
fprintf(wfile,'%s\n',regmodstr);
fprintf(wfile,'%s\n',lambdastr);
fprintf(wfile,'%s\n','3 INFMOD');
fprintf(wfile,'%s\n','-10. LAMBST');
fprintf(wfile,'%s\n','1. LAMBSP');
fprintf(wfile,'%s\n','12. LAMBRA');
fprintf(wfile,'%s\n','1.e-8 LAMBPR');
fprintf(wfile,'%s\n','200 LAMBIT');
fprintf(wfile,'%6.1f%s\n',15.0,' EXCI');
fclose(wfile);

