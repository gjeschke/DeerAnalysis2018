function [r,distr,dlow,dhigh]=eval_ftikreg_r(handles)
%
% Evaluate output data of FTIKREG (distance domain)


r=[];
distr=[];
dlow=[];
dhigh=[];
respath=which('ftikreg_r.sol'); % solution (output) data file
dataexist=exist(respath);
if dataexist~=2, return; end;

result=load('ftikreg_r.sol');

r=result(:,1);
distr=result(:,2);
dlow=result(:,3)-result(:,4);
dhigh=result(:,3)+result(:,4);
r=r';
distr=distr';
dlow=dlow';
dhigh=dhigh';




