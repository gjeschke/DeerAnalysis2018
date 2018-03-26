function [kernel,r,t,U,sm,X,V,L] = get_bas_Tikh(dimension,dt,rmin,rmax,deriv)
% [kernel,r,t,U,sm,X,V,L] = get_bas_Tikh(dimension,dt,rmin,rmax,deriv)
%
% On-the-fly computation of a kernel and its singular-value decomposition
% for use with Tikhonov regularization routine
%
% dimension     kernel dimension
% dt            time increment, defaults to 0.008 microseconds
% rmin          minimum distance, defaults to 1.5 nm
% rmax          maximum distance, defaults to 10 nm
% deriv         derivative order used in regularization, defaults to 2
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
end;

if ~exist('deriv','var') || isempty(deriv),
    deriv = 2;
end;

w0=2*pi*ny0; % angular frequencies
t=0:dt:(dimension-1)*dt; % time axis in µs
n=length(t); % length of time axis
r=linspace(rmin,rmax,dimension); % distance axis in nm
m=length(r); % length of distance axis
wd=zeros(1,m); % vector dipolar frequencies as a function of distance
kernel=zeros(m,n); % data array for kernel
% fprintf(1,'%i data points in time domain\n',n);
% fprintf(1,'%i data points in distance domain\n',m);
for k=1:m
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   wd(k)=wdd/(2*pi); % current dipolar frequency
%    if mod(k,10)==0,
%         fprintf(1,'%5.1f%% of fit table done.\n',100*k/m);
%    end;
   for l=0:1000 % 1000 averages over cos(theta) angle (powder average)
      x=l/1000; % current theta angle
      ww=wdd*(3*x^2-1); % dipolar frequency at current theta angle
      kernel(k,:)=kernel(k,:)+cos(ww*t); % add kernel contribution
   end
end

for k=1:m % loop form kernel normalization
   kernel(k,:)=kernel(k,:)/kernel(k,1); % normalize dipolar time evolution traces
end
kernel = kernel';
L = get_l(length(r),deriv); % differential operator matrix for derivative
[U,sm,X,V] = cgsvd(kernel,L);
