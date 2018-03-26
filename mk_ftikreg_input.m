function mk_ftikreg_input(handles)
%
% Convert experimental data format (deer_fit input) to internal formats
%
% fname     file name, output of hermit interpolation

if handles.non_negativity,
    regmodstr='222 REGMOD'; % was 223 before May 17th 2005
else
    regmodstr='122 REGMOD'; 
end;
lambdastr='1.0 LAMBDA'; % was 0.2 before May 17th 2005
if handles.regpar_input,
    if handles.non_negativity,
        regmodstr='212 REGMOD'; % was 213 before May 17th 2005
    else
        regmodstr='112 REGMOD'; 
    end;
    lambdastr=[num2str(handles.regpar) ' LAMBDA'];
end;
if handles.noise_averaging,
    regmodstr='112 REGMOD'; 
    lambdastr=[num2str(handles.regpar) ' LAMBDA'];
end;

% Range of analysis
%
rmin=handles.rmin; 
rmax=handles.rmax;
tdip=handles.A_tdip;
tmax=max(tdip);
dt=tdip(2)-tdip(1);
dipevo=handles.A_dipevo;
dipevo_err=handles.fit_rms_value*handles.A_dipevo_err*length(handles.A_dipevo_err)/sum(handles.A_dipevo_err);

% figure(11); clf;
% plot(tdip,dipevo_err,'k');

dospath=which('ftikreg_r_new.exe');
dospath=dospath(1:length(dospath)-length('ftikreg_r_new.exe'));

% New code Aug-7-2007 for logarithmic Tikhonov
% theoretically nice but numerically problematic
% if handles.log_Tikh,
%     logcluster=log(handles.A_cluster);
%     logdipevo=logcluster-log(handles.A_depth)*ones(size(logcluster));
% 
%     dipevo=logcluster;
%     figure(12); clf;
%     plot(tdip,dipevo,'k');
%     dospath=which('ftikreg_r_log.exe');
%     dospath=dospath(1:length(dospath)-length('ftikreg_r_log.exe'));
% end;

data=[tdip' dipevo'];
%data=[tdip' dipevo' dipevo_err'];
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
fprintf(wfile,'%6.1f%s\n',handles.bandwidth,' EXCI');
fclose(wfile);

