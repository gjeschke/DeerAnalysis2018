function [sim,sc]=deer_sim(r,distr,t,dipevo,exci)
% Simulation of DEER signal sim(t) from distance distribution distr(r)
% and scaling factor sc of distribution so that it fits original data set dipevo(t)
%
% (c) G. Jeschke, 2005
%
% r in nm
% t in microseconds
%

ny0=52.04; % dipolar frequency at 1 nm
level=1e-3*max(distr); % minimimum level for simulation of a point
sim=zeros(size(t));
rej=0;
mdepth=0;
fdepth=0;
for k=1:length(r),
    if distr(k)>level, % simulate only for distances with significant contribution
         nydd=ny0/r(k)^3; % dipolar frequency at current distance
         for m=1:1000, % theta loop for orientation dependence
             x=m/1000; % cos(theta) 
             nyac=nydd*(3*x^2-1); % dipolar frequency at current orientation
             warg=nyac/exci;
             weight=exp(-warg^2);
             fdepth=fdepth+distr(k);
             mdepth=mdepth+weight*distr(k);
             sim=sim+distr(k)*cos(2*pi*nyac*t)*weight; % add contribution to deer signal
         end;
    else
         rej=rej+1;
     end;
end;
sc=mdepth/fdepth;
sim=0.01*sim/mdepth;
sim=sim+0.99*ones(size(sim));

