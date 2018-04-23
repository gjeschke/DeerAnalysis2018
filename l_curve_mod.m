function [idx_corner,idx_AIC,rho,eta,reg_param] = l_curve_mod(K,L,S,noise)
%L_CURVE_MOD Solve Tikhonov for a range of regularization parameters and
% determine the optimal ones based on several criteria.

% Set defaults
if ~exist('noise','var'), noise = 0; end
unconstrainedTikhonov = true;

% Preparations
L = full(L);
nt = numel(S);
KtK = K.'*K;
LtL = L.'*L;
KtS = K.'*S;

% Get alpha range
reg_param = get_regparamrange(K,L,noise);

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
