function handles = compute_deernet(handles)
%
% Frame function that calls Kuprov/Worswick's deernet computation
% and stores the output in internal DeerAnalysis format
%
% G. Jeschke, 23.3.2018

set(handles.status_line,'String',sprintf('Computing DEERNet analysis [%s]...',handles.net_set));
set(gcf,'Pointer','watch');
drawnow;
silent = true;
netset_params = which(fullfile(handles.net_set,'netset_params.m'));
net_dir = fileparts(netset_params);

[~,ztpoi]=min(abs(handles.texp));

ttemp = handles.texp - handles.cutoff*ones(size(handles.texp))/1000;
[~,pcutoff]=min(abs(ttemp));

time_axis = handles.texp(ztpoi:pcutoff).';
deer_trace = real(handles.vexp(ztpoi:pcutoff)).';


[dist_axis,outcomes] = deernet(deer_trace,time_axis,net_dir,silent);
rexp = dist_axis.'/1000;
distr_ensemble = outcomes.';
[nm,~] = size(distr_ensemble);
distr = mean(distr_ensemble);
sc = 1/sum(distr);
distr = sc*distr;
handles.mean_distr = distr;
handles.distr_std = std(sc*distr_ensemble);
handles.A_r = rexp;
handles.A_distr = distr;
handles.A_low = min(sc*distr_ensemble);
handles.A_high= max(sc*distr_ensemble);

handles.A_deernet_t = time_axis.';
handles.A_tdip =  time_axis.';
handles.A_deernet_vexp = deer_trace.';
[handles.A_deernet_sim,handles.A_deernet_ff,handles.A_deernet_bckg] = fit_deernet_primary(handles,rexp,distr,time_axis.',deer_trace.');
handles.A_depth = 1 - handles.A_deernet_bckg(1);
handles.updated = 1;

deernet_ensemble_sim = zeros(nm,length(handles.A_deernet_sim));
deernet_ensemble_ff = zeros(nm,length(handles.A_deernet_ff));
deernet_ensemble_bckg = zeros(nm,length(handles.A_deernet_bckg));

for km = 1:nm
    [deernet_ensemble_sim(km,:),deernet_ensemble_ff(km,:),deernet_ensemble_bckg(km,:)] = fit_deernet_primary(handles,rexp,distr_ensemble(km,:),time_axis.',deer_trace.');
end

handles.deernet_ensemble_sim = deernet_ensemble_sim;
handles.deernet_ensemble_ff = deernet_ensemble_ff;
handles.deernet_ensemble_bckg = deernet_ensemble_bckg;

bmode = get_bckg_mode(handles);
if ~strcmp(bmode,'d')
    handles.old_bckg_mode = bmode;
end

set(handles.radiobutton_deernet_bckg,'Value',1);
set(handles.checkbox_deernet_error_bckg,'Value',1);
set(handles.radiobutton_deernet_bckg,'Enable','on');
set(handles.checkbox_deernet_error_bckg,'Enable','on');


handles.bckg_request_h = 0;
handles.bckg_request_e = 0;
handles.bckg_request_p = 0;
handles.bckg_request_d = 1;
set_bckg_mode(handles,'d');
handles.comp_net_set = handles.net_set;
handles.new_bckg = 1;
handles.new_distr = 1;

set(handles.status_line,'String','DEERNet analysis finished.');
set(gcf,'Pointer','arrow');
drawnow


