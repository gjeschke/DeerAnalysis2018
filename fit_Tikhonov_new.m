function handles=fit_Tikhonov_new(handles)
%
% Tikhonov regularization, with or without excitation bandwidth correction
% and with or without computation of a whole L curve
% based on Regularization Tools by Per Christian Hansen for unconstrained
% regularization and non-negative least squares fitting for constrained
% regularization
%
% G. Jeschke, 2015
%
% Modified by H.C. Hyde, 2011 
%  *Added feature to get number of Tikhonov free parameters and store in handles



LFlag=get(handles.select_L_curve,'Value');
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
if LFlag,
    msg=sprintf('Computing Tikhonov regularization L curve and result at corner'); 
else
    msg=sprintf('Computing Tikhonov regularization at alpha = %12.4f',handles.regpar);
end
set(handles.status_line,'String',msg);
set(handles.main_figure,'Pointer','watch');
drawnow;
if LFlag,
    [r,distr,rho,eta,reg_param,corner] = get_Tikhonov_new(handles);
else
    [r,distr,rho,eta,reg_param,corner] = get_Tikhonov_new(handles,handles.regpar);
end;
if LFlag,
	handles.regpar_vector = reg_param; % store Lcurve 
    handles.Lcurve_rho = rho;
    handles.Lcurve_eta = eta;
    handles.regpar = reg_param(corner);
else
    handles.regpar = reg_param;
end;
handles.A_r = r;
pstr=num2str(handles.regpar);
set(handles.regpar_edit,'String',pstr);
set(handles.status_line,'String','Simulating DEER data...');
if exflag || length(handles.A_tdip) > 1024,
    [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.A_cluster,handles.bandwidth);
else
    sim=get_td_fit(handles,r,distr);
    sc = 1;
end;
handles.moddepth_suppression=sc;
handles.A_sim=sim;
handles.A_distr=distr';
handles.A_low = distr';
handles.A_high = distr';
if length(r)*length(distr)>0,
    handles.updated=1;
end;
set(handles.status_line,'String','Distance domain Tikhonov regularization succeeded.');
modsim=ones(size(sim))-sim;
modexp=ones(size(handles.A_cluster))-handles.A_cluster;
sc=sum(modexp.*modexp)/sum(modsim.*modexp);
sim1=ones(size(modsim))-sc*modsim;
axes(handles.dipolar_evolution);
cla;
plot(handles.A_tdip,sim1,'r','LineWidth',1.5);
plot(handles.A_tdip,handles.A_cluster,'k');
pstr=num2str(handles.regpar);
set(handles.regpar_edit,'String',pstr);
set(handles.main_figure,'Pointer','arrow');
drawnow
diff0=handles.A_cluster-sim1;
deriv2=diff(distr,2);
if LFlag,
    set(handles.L_curve,'Enable','on');
    set(handles.L_curve,'Value',1);
    handles.Lcurve_distr=handles.A_distr;
    handles.Lcurve_sim=handles.A_sim;
    handles.regpars=handles.regpar_vector;
    handles.regpar_sel=corner;
    handles.regpar_opt=corner;
else
    set(handles.L_curve,'Enable','off');
    set(handles.L_curve,'Value',0);
    handles.Lcurve_distr=handles.A_distr;
    handles.Lcurve_sim=handles.A_sim;
    handles.regpars=handles.regpar;
    handles.regpar_sel=1;
    handles.regpar_opt=1;
end;
handles.mask=ones(size(handles.A_distr));
pstr=num2str(handles.regpar);
set(handles.regpar_edit,'String',pstr);
handles.updated=1;
set(handles.validate_Tikhonov,'Enable','on');
