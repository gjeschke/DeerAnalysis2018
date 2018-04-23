function [idx_corner,idx_AIC,rho,eta,reg_param] = l_curve_mod(K,L,S,noise)
%L_CURVE_MOD Solve Tikhonov for a range of regularization parameters and
% determine the optimal ones based on several criteria.

% Set defaults
%-------------------------------------------------------------
if ~exist('noise','var'), noise = 0; end
lgregpar_inc = 0.1;  % resolution of log10(alpha)
unconstrainedTikhonov = true;

L = full(L);
nt = numel(S);
KtK = K.'*K;
LtL = L.'*L;
KtS = K.'*S;

% Set alpha range
%-------------------------------------------------------------
minmax_ratio = 16*eps*1e6;  % Max. ratio of smallest to largest alpha

% The following scaling of the alpha range improves L curve corner detection
% for DEER applications.
minmax_ratio = minmax_ratio*2^(noise/0.0025);

% Get singular values
sv = gsvd(K,L,0);
sv = sv(end-2:-1:1); % sort in decreasing order
lgregpar_max = log10(sv(1));
lgregpar_min = log10(max([sv(p),sv(1)*minmax_ratio]));
lgregpar_max = floor(lgregpar_max/lgregpar_inc)*lgregpar_inc;
lgregpar_min = ceil(lgregpar_min/lgregpar_inc)*lgregpar_inc;
lgregpar = lgregpar_max:-lgregpar_inc:lgregpar_min;
reg_param = 10.^lgregpar;

% Calculate all metrics over alpha range
%-------------------------------------------------------------
for a = numel(reg_param):-1:1
  Q = KtK + reg_param(a)^2*LtL;
  if unconstrainedTikhonov
    P = Q\KtS; % unconstrained Tikhonov solution
  else
    P = fnnls(Q,KtS); % non-negative Tikhonov solution
  end
  Serr = K*P - S; % time-domain fit residuals
  H = K*(Q\K.'); % influence/projection/hat matrix
  rho(a) = norm(Serr); % time-domain fit error
  eta(a) = norm(L*P); % distance-domain roughness
  AICmetric(a) = nt*log(norm(Serr)^2/nt) + 2*trace(H);
  GCVmetric(a) = norm(Serr)^2/(1-trace(H)/nt)^2;
end

% Find alpha values at AIC and GCV metric minima
%-------------------------------------------------------------
[~,idx_AIC] = min(AICmetric);
[~,idx_GCV] = min(GCVmetric);

% Determine corner of L-curve
%-------------------------------------------------------------
eta = log(eta);
rho = log(rho);
etamin = min(eta);
rhomin = min(rho);
etamax = max(eta);
rhomax = max(rho);
eta_scaled = (eta-etamin)/(etamax-etamin);
rho_scaled = (rho-rhomin)/(rhomax-rhomin);
Lcorner_metric = eta_scaled.^2+rho_scaled.^2;
[~,idx_corner] = min(Lcorner_metric);

% Print results
%-------------------------------------------------------------
printResults = true;
if printResults
  fprintf('Regularization parameter search:\n');
  fprintf('  range:    alpha = %g to %g  log10(alpha) = %g to %g\n',...
    reg_param(1),reg_param(end),log10(reg_param(1)),log10(reg_param(end)));
  fprintf('  AIC:      alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    reg_param(idx_AIC),log10(reg_param(idx_AIC)),idx_AIC);
  fprintf('  GCV:      alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    reg_param(idx_GCV),log10(reg_param(idx_GCV)),idx_GCV);
  fprintf('  L-corner: alpha = %g  log10(alpha) = %g  (idx %d)\n',...
    reg_param(idx_corner),log10(reg_param(idx_corner)),idx_corner);
end
