
%% Manufactured unsteady solution
syms x y t g alpha beta vmax real
syms vx(y) u0 uinfty(x,y) u(x,y,t) 
vx = vmax*y;
u0 = g;
uinfty = g*(exp(vx*x/alpha) - 1)/(exp(vx/alpha) - 1);
u = u0 + (uinfty - u0)*(1 - exp(-beta*t^2));
u = simplify(u);
display(u)
%%
u_y0 = limit(u, y, 0);
display(u_y0)
%% h
syms h_y0(x,t) h_x0(y,t) h_x1(y,t)

h_x0 = -alpha*diff(u, x);
h_x0 = subs(h_x0, x, 0);
display(h_x0);

h_y0 = -alpha*diff(u, y);
h_y0 = limit(h_y0, y, 0);
display(h_y0);

h_y1 = alpha*diff(u, y);
h_y1 = subs(h_y1, y, 1);
display(h_y1);
%% s
syms s(x,y,t)
s = diff(u, t) + vx*diff(u, x) - ...
    alpha*(diff(diff(u, x), x) + diff(diff(u, y), y));
display(s)
s_y0 = limit(s, y, 0);
display(s_y0)