function [r,distr,sim,rho,eta,Tikhonov_logfile]=get_Tikhonov(t,ff,regpar)

mysystem=computer; % determine type of computer for selecting the proper Tikhonov executable
if strcmp(mysystem,'PCWIN') || strcmp(mysystem,'PCWIN64')
    win_flag=1;
else
    win_flag=0;
end;
mac_flag=ismac; % find out if the computer is a Mac
if win_flag % set proper filename for executable
    Tikh_call='ftikreg_r_old.exe';
else
    if mac_flag,
        Tikh_call='ftikreg_r_old.maci';
    else
        Tikh_call='ftikreg_r_old.out';
    end;
end;
mk_ftikreg_input_std_alone(t,ff,regpar);
callpath=which(Tikh_call);
dospath=callpath(1:length(callpath)-length(Tikh_call));
currdir=pwd;
fprintf(1,'%s%g%s\n','Alpha= ',regpar,'. Please wait while FTIKREG is executing...'); 
cd(dospath);
if win_flag
    [status,result]=dos(callpath);
else
    [status,result]=unix(callpath);
end;
Tikhonov_logfile=result;
save Tikhonov_log status result
cd(currdir);
if ~isempty(findstr('FTIKREG ENDED',result)),
    [r,distr]=eval_ftikreg_r_std_alone;
    fprintf(1,'Simulating DEER data...\n');
    sim=get_td_fit_std_alone(r,distr,t);
    if length(r)*length(distr)<=0,
        r=[];
        distr=[];
    end;
    diff0=ff-sim;
    deriv2=diff(distr,2);
    rho=log(sum(diff0.*diff0));
    eta=log(sum(deriv2.*deriv2));
else
    fprintf(1,'Distance domain Tikhonov regularization failed.\n');
    r=[];
    distr=[];
    sim=[];
end;

function [r,distr,dlow,dhigh]=eval_ftikreg_r_std_alone
%
% Evaluate output data of FTIKREG (distance domain)


r=[];
distr=[];
dlow=[];
dhigh=[];
dataexist=exist('ftikreg_r.sol','file'); % solution (output) data file
if ~dataexist, return; end;

result=load('ftikreg_r.sol');

r=result(:,1);
distr=result(:,2);
dlow=result(:,3)-result(:,4);
dhigh=result(:,3)+result(:,4);
r=r';
distr=distr';
dlow=dlow';
dhigh=dhigh';

