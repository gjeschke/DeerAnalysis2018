function void=make_bas();
%
% Computes kernel data for the DEER fit program
%
% G. Jeschke, 2002
%
ny0=52.04; % dipolar frequency at 1 nm for g=ge
w0=2*pi*ny0; % angular frequencies
t=0:0.008:8.192-0.008; % time axis in µs
n=length(t); % length of time axis
r=1.5:0.02:40; % distance axis in nm
m=length(r); % length of distance axis
wd=zeros(1,m); % vector dipolar frequencies as a function of distance
base=zeros(m,n); % data array for kernel
tic, % initialize computation time clock
for k=1:m,
   rk=r(k); % current distance
   wdd=w0/rk^3; % current dipolar angular frequency
   wd(k)=wdd/(2*pi); % current dipolar frequency
   pr=sprintf('%5.1f%s',100*k/m,'% of fit table done.');
   disp(pr); % display progress made in computation
   for l=1:1000, % 1000 averages over theta angle (powder average)
      th=l*pi/2000; % current theta angle
      st=sin(th);
      ct=cos(th);
      ww=wdd*(3*ct^2-1); % dipolar frequency at current theta angle
      base(k,:)=base(k,:)+st*cos(ww*t); % add kernel contribution
   end;
end;
toc, % output of computation time required
for k=1:m, % loop form kernel normalization
   base(k,:)=base(k,:)/base(k,1); % normalize dipolar time evolution traces
end;
save('pake_base40.mat','base','r','t','wd'); % save kernel data
void=1;
