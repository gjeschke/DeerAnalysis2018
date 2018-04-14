function [K,r,t] = get_bas_Tikh(dimension,dt,rmin,rmax)
% [K,r,t] = get_bas_Tikh(dimension,dt,rmin,rmax)
%
% Compute dipolar kernel for use with Tikhonov regularization
%
% Input:
%   dimension     kernel dimension
%   dt            time increment, defaults to 0.008 microseconds
%   rmin          minimum distance, defaults to 1.5 nm
%   rmax          maximum distance, defaults to 10 nm
%
% Output:
%   r             distance vector, nm
%   t             time vector, microseconds
%   K             kernel

ny0=52.04; % dipolar frequency at 1 nm for g=ge, MHz

if ~exist('dt','var') || isempty(dt)
    dt = 0.008;
end

if ~exist('rmin','var') || isempty(rmin)
    % minimum detectable distance, the largest dipolar frequency (shoulders of
    % the Pake doublet is 85% of the Nyquist frequency)
    rmin = (4*dt*ny0/0.85)^(1/3); 
end

if ~exist('rmax','var') || isempty(rmax)
    % maximum detectable distance, 2 microseconds dipolar evolution time
    % corresponds to 6 nm
    rmax = 6*(dimension*dt/2)^(1/3);
end

t = (0:dimension-1).'*dt; % time axis in µs
nt = length(t); % length of time axis
r = linspace(rmin,rmax,dimension); % distance axis in nm
nr = length(r); % length of distance axis
K = zeros(nt,nr); % kernel matrix
wdd = 2*pi*ny0./r.^3; % perpendicular dipolar angular frequencyies for all r

Method = 2;
switch Method
  case 1
    % Method using explicit numerical powder average (slow)
    %----------------------------------------------------------
    ntheta = 1001;
    costheta = linspace(0,1,ntheta);
    for ir = 1:nr
      Kr = 0;
      for l = 1:ntheta % average over cos(theta) angle (powder average)
        ww = wdd(ir)*(1-3*costheta(l)^2); % dipolar frequency at current theta
        Kr = Kr + cos(ww*t); % add kernel contribution
      end
      % normalize dipolar time evolution trace
      K(:,ir) = Kr/Kr(1);
    end
  case 2
    % Method using Fresnel integrals (fast)
    %----------------------------------------------------------
    % using John D'Errico's implementation of the Fresnel integrals
    % from the Matlab file exchange
    % Equations: see Edwards/Stoll, J.Magn.Reson. 270, 87-97 (2016), Eq.(6)
    % https://doi.org/10.1016/j.jmr.2016.06.021
    for ir = 1:nr
      ph = wdd(ir)*abs(t);
      y = sqrt(6*ph/pi);
      K(:,ir) = (cos(ph).*fresnelC(y)+sin(ph).*fresnelS(y))./y; % div by zero for y=0
    end
    % Correct at t = 0, since the Fresnel expression
    % involves division by zero for t = 0.
    K(t==0,:) = 1;
  otherwise
    error('Unknown calculation method for dipolar kernel.')
end
