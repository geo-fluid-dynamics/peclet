syms x t g v alpha beta real
syms ue(x) ue_v0(x) um(x,t) um_v0(x) sm(x,t) sm_v0(x,t) hm(t) hm_v0(t)
%% Exact steady solution
ue = g*(exp(v*x/alpha) - 1)/(exp(v/alpha) - 1);
ue_v0 = limit(ue, v, 0);
display(ue)
display(ue_v0)
%% Manufactured unsteady solution
um = g*(1 + ((exp(v*x/alpha) - 1)/(exp(v/alpha) - 1) - 1)*...
    (1 - exp(-beta*t^2)));
um_v0 = limit(um, v, 0);
display(um)
display(um_v0)
%% Manufactured Neumann boundary conditions
umx = diff(um, x);
hm = subs(-alpha*umx, x, 0);
display(hm)
hm_v0 = limit(hm, v, 0);
display(hm_v0)
%% Manufactured source term
umt = diff(um, t);
umxx = diff(umx, x);
sm = umt + v*umx - alpha*umxx;
display(sm)
sm_v0 = limit(sm, v, 0);
display(sm_v0)
