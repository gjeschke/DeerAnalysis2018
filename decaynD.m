function sim=decaynD(v,x,hom_dim,handles)
%
%
% distr=handles.Pake_r.^(handles.hom_dim-1);
% distr=0.01*distr/sum(distr);
% logB=distr*handles.Pake_kernel;
fac=1;
if max(x)>max(handles.Pake_t),
    fac=max(x)/max(handles.Pake_t);
end;

% sim0=v(2)*exp(fac*16000*v(1)*logB);
% sim=interp1(handles.Pake_t,sim0,x/fac,'pchip');

sim=real(v(2)*exp(-(v(1)*x).^(hom_dim/3)));
% keyboard
% disp('Hallo');
% sim=v(2)*exp(-v(1)*x.^(n/3));
