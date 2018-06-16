function handles=fit_Tikhonov(handles)
%
% Tikhonov regularization, with or without excitation bandwidth correction
% and with or without computation of a whole L curve
% calls slightly modified versions of the FORTRAN program FTIKREG
%
% G. Jeschke, 2006
%
% Modified by H.C. Hyde, 2011 
%  *Added feature to get number of Tikhonov free parameters and store in handles


mysystem=computer; % determine type of computer for selecting the proper Tikhonov executable
if strcmp(mysystem,'PCWIN') || strcmp(mysystem,'PCWIN64')
    win_flag=1;
else
    win_flag=0;
end;
mac_flag=ismac; % find out if the computer is a Mac
exflag=get(handles.exci_bandwidth_corr,'Value'); % check, if excitation bandwidth correction is selected
if win_flag % set proper filename for executable
    Tikh_call='ftikreg_r_old.exe';
else
    if mac_flag,
        Tikh_call='ftikreg_r_old.maci';
    else
        Tikh_call='ftikreg_r_old.out';
    end;
end;
if exflag, % change executable name, if excitation bandwidth correction was selected
    Tikh_call='ftikreg_r_new.out';
    if win_flag
        Tikh_call='ftikreg_r_new.exe';
    end;
    if mac_flag,
        Tikh_call='ftikreg_r_new.maci';
    end;        
end;
Lflag=get(handles.select_L_curve,'Value'); % check, if whole L curve should be computed
if Lflag,
	regpar_vector=handles.rpv; % for L curve computation, a whole vector of regulation parameters exist
else
    regpar_vector=handles.regpar; % otherwise, the "vector" is a single value
end;
nn=length(regpar_vector);
regdistr=zeros(nn,4096);
reglow=regdistr;
reghigh=regdistr;
regsim=zeros(nn,length(handles.A_tdip));
rho=zeros(size(regpar_vector));
eta=zeros(size(regpar_vector));
handles.comp_r_min=handles.rmin;
handles.comp_r_max=handles.rmax;
status_figure('L curve: Close to stop.');
for rr=1:length(regpar_vector),
    handles.regpar=regpar_vector(rr);
    handles.regpar_input=1;
	mk_ftikreg_input(handles);
	callpath=which(Tikh_call);
	dospath=callpath(1:length(callpath)-length(Tikh_call));
	currdir=pwd;
    msg=sprintf('%s%g%s','Alpha= ',handles.regpar,'. Please wait while FTIKREG is executing...'); 
	set(handles.status_line,'String',msg);
	set(handles.main_figure,'Pointer','watch');
    drawnow;
	cd(dospath);
    if exist('ftikreg_r.sol','file'),
        delete('ftikreg_r.sol');
    end;

    if win_flag
        [status,result]=dos(callpath);
    else
        [status,result]=unix(callpath);
    end;
    handles.Tikhonov_logfile=result;
    
    % New code (HCH): Get "NUMBER OF FREE PARAMETERS" from Tikhonov fit results log
    % and interpret this value as the effective number of free parameters. It
    % is equivalent to the number of 'r' values at which P(r) is defined [Ref. 1].
    % 1) G. Jeschke, et al. DeerAnalysis2006 - A comprehensive software package
    %    for analyzing pulsed ELDOR data. Appl. Magn. Reson. 30, 473–498 (2006).
    % --------
    Istr = strfind(result,'CAFSQP > NUMBER OF FREE PARAMETERS');  %first index of target line in result log
    if ~isempty(Istr)
         handles.Tikh_nfp_vector(rr) = round(str2num(result(Istr+38:Istr+42)));
    else handles.Tikh_nfp_vector = []; 
    end
    % --------

    save Tikhonov_log status result
	cd(currdir);
	set(handles.main_figure,'Pointer','arrow');
    drawnow;
	path=handles.project_dir;
	bas_name=handles.bas_name;
	if length(findstr('FTIKREG ENDED',result))>0,
%         disp('Before');
%         disp(handles.regpar);
        [r,distr,dlow,dhigh]=eval_ftikreg_r(handles);
        handles.A_r=r;
        sc=1;
        axes(handles.distance_distribution);
        set(handles.distance_distribution,'NextPlot','replace');
	    cla;
%         disp('0');
%         disp(handles.regpar);
        set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
        plot(r,distr,'k');       
		set(handles.status_line,'String','Simulating DEER data...');
		set(handles.main_figure,'Pointer','watch');
        if exflag,
            [sim,sc]=deer_sim(r,distr,handles.A_tdip,handles.A_cluster,handles.bandwidth);
        else,
            sim=get_td_fit(handles,r,distr);
        end;
        handles.moddepth_suppression=sc;
        regsim(rr,:)=sim;
        handles.A_sim=sim;
        handles.A_distr=distr;
        handles.A_low=dlow;
        handles.A_high=dhigh;
        regdistr(rr,1:length(distr))=distr;
        reglow(rr,1:length(distr))=dlow;
        reghigh(rr,1:length(distr))=dhigh;
%        disp(result);
        poivec=findstr('REGULARIZATION PARAMETER    :',result);
        if ~isempty(poivec),
            poi=poivec+length('REGULARIZATION PARAMETER    :');
            answer=result(poi:length(result));
            answer=strtok(answer);
            handles.regpar=str2num(answer);
        end;
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
        set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
		set(handles.main_figure,'Pointer','arrow');
        drawnow
        diff0=handles.A_cluster-sim1;
        deriv2=diff(distr,2);
		rho(rr)=log(sum(diff0.*diff0));
        eta(rr)=log(sum(deriv2.*deriv2));
    else
        set(handles.status_line,'String','Distance domain Tikhonov regularization failed.');
        set(handles.regpar_edit,'String','n.a.');
        handles.updated=0;
        handles.A_r=[];
        handles.A_distr=[];
        handles.A_sim=[];
	end;
    comp_status=status_figure(rr/length(regpar_vector));     
    if ~comp_status, regpar_vector=regpar_vector(1:rr); break; end;
end;
drawnow;
if comp_status, status_figure(1); end;
% disp('2');
% disp(handles.regpar);
if length(regpar_vector)>1, % compute L curve and its corner
	regdistr=regdistr(:,1:length(handles.A_r));
	reglow=reglow(:,1:length(handles.A_r));
	reghigh=reghigh(:,1:length(handles.A_r));
	handles.Lcurve_rho=rho;
	handles.Lcurve_eta=eta;
	poi=get_l_corner(rho,eta);
	handles.regpar_sel=poi;
    handles.regpar_opt_Lc=poi;
    handles.regpar=regpar_vector(poi);
    handles.regpars=regpar_vector;
	opt_reg=regpar_vector(poi);
	handles.A_distr=regdistr(poi,:);
    handles.A_low=reglow(poi,:);
    handles.A_high=reghigh(poi,:);
	handles.A_sim=regsim(poi,:);
	handles.regpar_input=1;
	handles.regpar=opt_reg;
    handles.Lcurve_distr=regdistr;
    handles.Lcurve_low=reglow;
    handles.Lcurve_high=reghigh;
    handles.Lcurve_sim=regsim;
    set(handles.L_curve,'Enable','on');
    set(handles.L_curve,'Value',1);
else
    set(handles.L_curve,'Enable','off');
    set(handles.L_curve,'Value',0);
    handles.Lcurve_distr=handles.A_distr;
    handles.Lcurve_sim=handles.A_sim;
    handles.regpars=handles.regpar;
    handles.regpar_sel=1;
    handles.regpar_opt_Lc=1;
end;
handles.mask=ones(size(handles.A_distr));
% disp('3');
% disp(handles.regpar);
set(handles.regpar_edit,'String',num2str(handles.regpar,handles.regpar_edit_strformat));
handles.updated=1;
set(handles.validate_Tikhonov,'Enable','on');

