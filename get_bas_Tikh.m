function [kernel,r,t] = get_bas_Tikh(dimension,dt,rmin,rmax)
% [kernel,r,t] = get_bas_Tikh(dimension,dt,rmin,rmax)
%
% Computation of a kernel for use with Tikhonov regularization routine
%
% dimension     kernel dimension
% dt            time increment, defaults to 0.008 microseconds
% rmin          minimum distance, defaults to 1.5 nm
% rmax          maximum distance, defaults to 10 nm
%
% G. Jeschke, 2015
%

ny0=52.04; % dipolar frequency at 1 nm for g=ge

if ~exist('dt','var') || isempty(dt)
    dt = 0.008;
end

if ~exist('rmin','var') || isempty(rmin)
    rmin = (4*dt*ny0/0.85)^(1/3); % minimum detectable distance, 
                                  % the largest dipolar frequency
                                  % (shoulders of the Pake doublet is 85%
                                  % of the Nyquist frequency)
end

if ~exist('rmax','var') || isempty(rmax)
    rmax = 6*(dimension*dt/2)^(1/3); % maximum detectable distance, 
                                     % 2 microseconds dipolar evolution
                                     % time correspond to 6 nm
end

w0=2*pi*ny0; % angular frequencies
t=(0:dimension-1).'*dt; % time axis in µs
nt=length(t); % length of time axis
r=linspace(rmin,rmax,dimension); % distance axis in nm
nr=length(r); % length of distance axis
%wd=zeros(1,nr); % vector dipolar frequencies as a function of distance
kernel=zeros(nt,nr); % kernel matrix
ntheta = 1001;
costheta = linspace(0,1,ntheta);
for k=1:nr
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   %wd(k)=wdd/(2*pi); % current dipolar frequency
   for l=1:ntheta % average over cos(theta) angle (powder average)
      ww=wdd*(3*costheta(l)^2-1); % dipolar frequency at current theta angle
      kernel(:,k)=kernel(:,k)+cos(ww*t); % add kernel contribution
   end
   % normalize dipolar time evolution trace
   kernel(:,k) = kernel(:,k)/kernel(1,k);
end
