function [sim,sc]=deer_sim(r,distr,t,dipevo,exci)
% Simulation of DEER signal sim(t) from distance distribution distr(r)
% and scaling factor sc of distribution so that it fits original data set dipevo(t)
%
% (c) G. Jeschke, 2005
%
% r in nm
% t in microseconds
%

ny0=52.04; % dipolar frequency at 1 nm (MHz)
level=1e-3*max(distr); % minimimum level for simulation of a point
sim=zeros(size(t));
mdepth=0;
fdepth=0;
neglect = distr<=level;
for k=1:length(r)
    % simulate only for distances with significant contribution
    if neglect(k), continue; end 
    
    nydd=ny0/r(k)^3; % dipolar frequency at current distance
    for m=1:1000 % theta loop for orientation dependence
      costheta=m/1000; % cos(theta)
      nyac=nydd*(3*costheta^2-1); % dipolar frequency at current orientation
      weight=exp(-(nyac/exci)^2);
      mdepth=mdepth+weight*distr(k);
      fdepth=fdepth+distr(k);
      sim=sim+distr(k)*cos(2*pi*nyac*t)*weight; % add contribution to deer signal
    end
end
sim=0.01*sim/mdepth;
sim=sim+0.99*ones(size(sim));
sc=mdepth/fdepth;

