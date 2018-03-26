function [sim,ff,bckg] = fit_deernet_primary(handles,rexp,distr,texp,vexp)

simff = get_td_fit(handles,rexp,distr);

dipevo = 100*(simff-0.99);

depth0 = handles.A_depth;
k0 = -log(vexp(end)/(1-depth0))/texp(end);
v0 = [3 k0 depth0];

v = fminsearch(@rms_net_bckg,v0,[],texp,vexp,dipevo);
[sim,ff,bckg] = sim_net_bckg(v,texp,dipevo);


function [sim,ff,bckg] = sim_net_bckg(v,texp,simff)

hom_dim = v(1);
kdec = v(2);
moddepth = v(3);

bckg = real(exp(-(kdec*texp).^(hom_dim/3)));
ff = moddepth*simff + (1-moddepth);
sim = ff.*bckg;
bckg = bckg*(1-moddepth);

function rmsd = rms_net_bckg(v,texp,vexp,simff)

sim = sim_net_bckg(v,texp,simff);
rmsd = sqrt(sum((vexp-sim).^2));
