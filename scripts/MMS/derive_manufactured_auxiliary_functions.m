syms x t g a Pe_r beta real
syms ue(x) ue_a0(x) um(x,t) um_a0(x) sm(x,t) sm_a0(x,t) hm(t) hm_a0(t)
%% Exact steady solution
ue = g*(exp(Pe_r*a*x) - 1)/(exp(Pe_r*a) - 1);
ue_a0 = limit(ue, a, 0);
display(ue)
display(ue_a0)
%% Manufactured unsteady solution
um = g*(1 + ((exp(Pe_r*a*x) - 1)/(exp(Pe_r*a) - 1) - 1)*...
    (1 - exp(-beta*t^2)));
um_a0 = limit(um, a, 0);
display(um)
display(um_a0)
%% Manufactured Neumann boundary conditions
umx = diff(um, x);
hm = subs(-1/Pe_r*umx, x, 0);
display(hm)
hm_a0 = limit(hm, a, 0);
display(hm_a0)
%% Manufactured source term
umt = diff(um, t);
umxx = diff(umx, x);
sm = umt + a*umx - 1/Pe_r*umxx;
display(sm)
sm_a0 = limit(sm, a, 0);
display(sm_a0)
