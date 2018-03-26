function [bckg,handback]=fit_bckg(handles,texp,t_fit,td_fit)
%
% Fits background decay according to the model selected by
% handles.background
% 0  fractal, n variable,  exp(-k*t^(n/3))
% 1  n-dimensional, n fixed, exp(-k*t^(n/3))
% 2  three-dimensional, exp(-k*t)
% 3  polynomial
% 4  user-defined function in handles.bckg_fct or numerical background in handles.bckg_data
% 5  homogeneous background 


man_bckg_flag=get(handles.manual_bckg,'Value');

pflag=get(handles.bckg_poly,'Value');
if pflag, back_model=3; end;
hflag=get(handles.bckg_homogeneous,'Value');
if hflag,
    dflag=get(handles.bckg_fit_dim,'Value');
    if dflag,
        back_model=0;
    else,
        back_model=1;
    end;
    if ~dflag && handles.hom_dim==3,
        back_model=2;
    end;
    if man_bckg_flag,
        back_model=5;
    end;
end;
uflag=get(handles.bckg_exp,'Value');
if uflag, 
    back_model=4;
    set(handles.density_text,'String','Rel. dens.');
else,
    set(handles.density_text,'String','Density');
end;    
        
[poly,s]=polyfit(t_fit,log(td_fit),1); % linear fit of logarithm
v0=[-poly(1) 1];
switch back_model,
    case 0,
        v0=[-poly(1) 1 3];
        v1=fminsearch(@rms_stretched,v0,[],t_fit,td_fit,handles);
        % bckg=decay_stretched(v1,texp);
%         distr=handles.Pake_r.^(v1(3)-1);
%         distr=0.01*distr/sum(distr);
%         logB=distr*handles.Pake_kernel;
        bckg=decaynD(v1(1:2),texp,v1(3),handles);
        dens=v1(1);
        fdim=v1(3);
        pstr=sprintf('%5.2f',fdim);
        set(handles.bckg_dim_edit,'String',pstr);
        handles.hom_dim=fdim;
        handles.bckg_dens=dens;
    case 1,
        distr=handles.Pake_r.^(handles.hom_dim-1);
        distr=0.01*distr/sum(distr);
        logB=distr*handles.Pake_kernel;
        v1=fminsearch(@rmsnD,v0,[],t_fit,td_fit,logB,handles);
        bckg=decaynD(v1,texp,handles.hom_dim,handles);
        dens=v1(1);
        handles.bckg_dens=dens;
    case 2,
		bckg=exp(polyval(poly,texp)); % background is exponential of that
        dens=v0(1);
        handles.bckg_dens=dens;
    case 3,
        [poly,s]=polyfit(t_fit,log(td_fit),handles.poly_order);
        bckg=exp(polyval(poly,texp));
        handles.polynomial=poly;
        dens=-1;
    case 4,
        bckg0=exp(polyval(handles.polynomial,t_fit));
        v0=[1 0.8];
        v1=fminsearch(@rms_ubckg,v0,[],td_fit,bckg0);
        bckg1=exp(polyval(handles.polynomial,texp));
        bckg=v1(2)*exp(v1(1)*log(bckg1));
        dens=v1(1)/handles.calib_density;
        handles.bckg_dens=dens;
    case 5,
        %targ=texp.^(handles.hom_dim/3);
        v1(1)=handles.man_k;
        v1(2)=1;
        bckg=decaynD(v1,texp,handles.hom_dim,handles);
        bckg=(1-handles.man_depth)*bckg;
        % bckg=(1-handles.man_depth)*exp(-dens*targ);
        dens=handles.man_k;
end;
% figure(13); clf;
% hold on;
% plot(t_fit,td_fit,'k');
% plot(t_fit,bckg0,'g');
% plot(texp,bckg,'r');
% figure(14);
% bckg=real(bckg);
pstr=sprintf('%6.3f',dens*handles.calib_density);
if dens>=0,
	set(handles.bckg_density,'String',pstr);
else
	set(handles.bckg_density,'String','n.a.');
end;
handles.density_value=dens*handles.calib_density;

handback=handles;

bckg=real(bckg);
