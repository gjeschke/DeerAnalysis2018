function [rout,distr,rho,eta,reg_param,corneridx] = get_Tikhonov_new(handles,reg_param)

% Based on regularization tools by Per Christian Hansen
% see: http://www2.compute.dtu.dk/~pcha/Regutools/index.html
% as well as: http://ch.mathworks.com/matlabcentral/fileexchange/52-regtools

persistent kdim kernel r t U sm X V L

% Determine whether L curve should be calculated or not
calcLcurve = ~exist('reg_param','var') || length(reg_param)~=1;

tdip = handles.A_tdip;

if length(tdip) > 2048 % use standard kernel for very long data sets
    kernel = handles.Tikh_kernel;
    r = handles.Tikh_r;
    t = handles.Tikh_t;
    U = handles.Tikh_U;
    sm = handles.Tikh_sm;
    X = handles.Tikh_X;
    V = handles.Tikh_V;
    L = handles.Tikh_L;
    kdim = length(r);
    n1 = length(t);
    tf = linspace(0,max(tdip),n1);
    dt2 = tf(2)-tf(1);
    ff = interp1(tdip,handles.A_dipevo,tf);
else
    if isempty(kdim) || kdim ~= length(tdip) % use previously stored kernel if dimension matches
        kdim = length(tdip);
        [kernel,r,t] = get_bas_Tikh(kdim);
        L = get_l(length(r),2); % differential regulariztion operator matrix
        [U,sm,X,V] = cgsvd(kernel,L);
    end
    ff = handles.A_dipevo;
    dt2 = tdip(2)-tdip(1);
end

dt = t(2)-t(1);
sc = (dt2/dt)^(1/3);

% % Code for testing
% figure(17); clf;
% plot(tdip,handles.A_dipevo,'k');
% hold on;
% plot(tf,ff,'r');
% fprintf(1,'New data size is (%i,%i).\n',size(ff));

if calcLcurve
    % Compute L curve rho and eta (without nonnegativity constraint) and its corner
    [corneridx,rho,eta,reg_param] = l_curve_mod(U,sm,ff','Tikh',L,V,handles.fit_rms_value);
    % distr0 = tikhonov(U,sm,X,ff',reg_param(corner));
% %   Code for testing
%     figure(8); clf;
%     plot(rho,eta,'k.');
%     hold on;
%     plot(rho(corner),eta(corner),'ro');
    calcAICmetric = false;
    if calcAICmetric
        unconstrainedAIC = true;
        % Calculate AIC metric for all alpha values in range
        S = ff(:);
        nt = numel(S);
        KtK = kernel.'*kernel;
        LtL = L.'*L;
        KtS = kernel.'*S;
        for a = numel(reg_param):-1:1
          alpha = reg_param(a);
          Q = KtK + alpha^2*LtL;
          if unconstrainedAIC
            P = Q\KtS; % unconstrained Tikhonov solution
          else
            P = fnnls(Q,KtS); % non-negative Tikhonov solution
          end
          Serr = kernel*P - S; % time-domain fit error
          H = kernel*(Q\kernel.'); % influence/projection/hat matrix
          AICmetric(a) = nt*log(norm(Serr)^2/nt) + 2*trace(H);
        end
        % Find alpha value at AIC metric minimum
        [~,idx_est] = min(AICmetric);
        fprintf('AIC optimum: alpha = %g  log10(alpha) = %g  (idx %d)\n',...
          reg_param(idx_est),log10(reg_param(idx_est)),idx_est);
    end
else
    % Compute rho and eta for single regularization parameter
    [~,rho,eta] = tikhonov(U,sm,X,ff',reg_param);
    corneridx = 1;
end

rout = sc*r;

% Solve non-negativity-constrained Tikhonov problem
%---------------------------------------------------------------
% solver 1: 2x matrix size, slow
%    argmin(||C*x-d||^2) with C = [K;alpha*L] and d = [ff;zeros]
% solver 2: 1x matrix size, fast
%    argmin(||ff-K*distr||^2+alpha^2*||L*P||^2), with functional multiplied out
nonnegSolver = 2;
alpha = reg_param(corneridx);

switch nonnegSolver
  case 1
    C = [kernel;alpha*L];
    % options = optimset('Display','off','TolX',10*eps*max(size(C))*norm(kernel,1));
    options = optimset('Display','off','TolX',1000*eps);
    [m,~] = size(L);
    d = [ff';zeros(m,1)];
    distr = lsqnonneg(C,d,options);
  case 2
    Q = (kernel.'*kernel) + alpha^2*(L.'*L);
    distr = fnnls(Q,kernel.'*ff(:));
end
