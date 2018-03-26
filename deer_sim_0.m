function sim =deer_sim_0(r,distr,t)
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
fdepth=0;
for k=1:length(r),
    if distr(k)>level, % simulate only for distances with significant contribution
         nydd=ny0/r(k)^3; % dipolar frequency at current distance
         for m=1:1000, % theta loop for orientation dependence
             x=m/1000; % cos(theta) 
             nyac=nydd*(3*x^2-1); % dipolar frequency at current orientation
             fdepth=fdepth+distr(k);
             sim=sim+distr(k)*cos(2*pi*nyac*t); % add contribution to deer signal
         end;
    else
         rej=rej+1;
     end;
end;
sim=0.01*sim/fdepth;
sim=sim+0.99*ones(size(sim));

