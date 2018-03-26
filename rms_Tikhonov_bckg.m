function rms=rms_Tikhonov_bckg(v,handles,texp,vexp,dim)
% computes the r.m.s.d. of the form factor computed from a Tikhonov-derived
% distance distribution from the experimental form factor
% a homogeneous background correction with dimension dim is applied, whose
% parameters are given in vector v
%
% v(1)  decay time constant
% v(2)  modulation depth
% texp  time axis of zero-time corrected, normalized and cutoff original data
% vexp  zero-time corrected and cutoff original data
%
% (c) G. Jeschke, 2008

targ=texp.^(dim/3);
bckg=(1-v(2))*exp(-v(1)*targ);
%         figure(13); clf;
%         plot(texp,vexp,'k');
%         hold on;
%         plot(texp,bckg,'r');
dipevo=real(vexp)-bckg;
dipevo=dipevo./bckg; % divide by background, eqn [13]
cluster=real(vexp)./bckg;
cluster=cluster/max(cluster);
handles.A_cluster=cluster;
handles.tdip=texp;
handles.dipevo=dipevo;
rms=get_Tikhonov(handles);
%         figure(14); clf;
%         plot(texp,cluster,'k');
%         hold on;
%         plot(texp,sim,'r');
pstr=sprintf('%s%6.3f%s%6.3f%s%6.4f','Optimizing bckg. fit parameters: k= ',v(1),', depth= ',v(2),' , r.m.s.d.: ',rms);
set(handles.status_line,'String',pstr);
% disp(pstr);
drawnow;
% keyboard

function rms=get_Tikhonov(handles,cluster),
%
% Tikhonov regularization, with or without excitation bandwidth correction
% and with or without computation of a whole L curve
% calls slightly modified versions of the FORTRAN program FTIKREG
%
% G. Jeschke, 2007
%


mysystem=computer;
if strcmp(mysystem(1:5),'PCWIN')
    win_flag=1;
else
    win_flag=0;
end;
sc=1;
exflag=get(handles.exci_bandwidth_corr,'Value');
if win_flag
    Tikh_call='ftikreg_r_old.exe';
else
    Tikh_call='ftikreg_r_old.out';
end;
if exflag,
    Tikh_call='ftikreg_r_new.out';
    if win_flag
        Tikh_call='ftikreg_r_new.exe';
    end;
end;
mk_ftikreg_input(handles);
callpath=which(Tikh_call);
dospath=callpath(1:length(callpath)-length(Tikh_call));
currdir=pwd;
msg=sprintf('%s%g%s','Alpha= ',handles.regpar,'. Please wait while FTIKREG is executing...'); 
% set(handles.status_line,'String',msg);
% if ~get(handles.statusbar_off,'Value'),
%     set(handles.Tikhonov_validation,'Pointer','watch');
%     drawnow;
% end;
cd(dospath);
if win_flag
    [status,result]=dos(callpath);
else
    [status,result]=unix(callpath);
end;
cd(currdir);
% if ~get(handles.statusbar_off,'Value'),
%     set(handles.Tikhonov_validation,'Pointer','arrow');
%     drawnow;
% end;
% path=main_handles.project_dir;
% bas_name=main_handles.bas_name;
if length(findstr('FTIKREG ENDED',result))>0,
    [r,distr,dlow,dhigh]=eval_ftikreg_r(handles);
    sig_distr=(dhigh-dlow)/2;
    exflag=get(handles.exci_bandwidth_corr,'Value');
    if exflag,
        [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.A_cluster,handles.bandwidth);
    else,
        sim=get_td_fit(handles,r,distr);
    end;
    modsim=ones(size(sim))-sim;
    modexp=ones(size(handles.A_cluster))-handles.A_cluster;
    sc=sum(modexp.*modexp)/sum(modsim.*modexp);
    sim1=ones(size(modsim))-sc*modsim;
    diff0=handles.A_cluster-sim1;
    rms=sqrt(sum(diff0.*diff0)/(length(diff0)-1));
else,
    rms=1e6;
end;