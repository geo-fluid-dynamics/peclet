syms x t g a Pe_r beta real
syms ue(x) ue_a0(x) um(x,t) um_a0(x) sm(x,t) sm_a0(x,t) hm(t) hm_a0(t)
%% Exact steady solution
ue = g*(1 + (a*Pe_r)/(1 - exp(a*Pe_r))/(a*Pe_r)*(exp(a*Pe_r) ...
    - exp(a*Pe_r*x)));
ue_a0 = limit(ue, a, 0);
display(ue)
display(ue_a0)
%% Manufactured unsteady solution
um = g*(1 + (exp(a*Pe_r) - exp(a*Pe_r*x))/(1 - exp(a*Pe_r))*...
    (1 - exp(-beta*t^2)));
um_a0 = limit(um, a, 0);
display(um)
display(um_a0)
%% Manufactured Neumann boundary conditions
umx = diff(um, x);
hm = subs(-umx, x, 0);
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
%% Plot for various a and Pe_r
aList = [0, -1, -5];
aColors = {'k', 'b', 'r'};
PeList = [1, 0.5];
PeStyles = {'-', '--'};
g = -1;
beta = 10;

uSteady_Figure = figure();
uSteady_FigureLegendStrings = {};

h_Figure = figure();
xlabel('$\tilde{t}$', 'Interpreter', 'latex');
        ylabel('$\tilde{h}$', 'Interpreter', 'latex');
h_FigureLegendStrings = {};

for ia = 1:length(aList)
    a = aList(ia);
    for iPe = 1:length(PeList)
        if a == 0 && iPe > 1
            continue
        end
        Pe_r = PeList(iPe);
        %% Make callable functions from symbolic functions
        if a == 0
            f_ue = matlabFunction(subs(ue_a0));
        else
            f_ue = matlabFunction(subs(ue));
        end

        if a == 0
            f_um = matlabFunction(subs(um_a0));
        else
            f_um = matlabFunction(subs(um));
        end

        if a == 0
            f_sm = matlabFunction(subs(sm_a0));
        else
            f_sm = matlabFunction(subs(sm));
        end

        if a == 0
            f_hm = matlabFunction(subs(hm_a0));
        else
            f_hm = matlabFunction(subs(hm));
        end
        %% Plot steady state solution vs exact steady state
        figure(uSteady_Figure)
        hold on
        xh = 0:0.01:1;
        plot(xh, f_ue(xh), ':k',...
            xh, f_um(1, xh), [PeStyles{iPe}, aColors{ia}])
        hold off
        legend({'Exact steady solution', 'Steady state of manufactured solution'})
        uSteady_FigureLegendStrings = [uSteady_FigureLegendStrings...
            {['Exact: $a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$']},...
            {['Manufactured: $a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$']}]; %#ok<AGROW>
        %% Plot exact solutions alone
        %% Plot steady state solution vs exact steady state
        figure(uSteady_Figure)
        hold on
        xh = 0:0.01:1;
        plot(xh, f_ue(xh), [PeStyles{iPe}, aColors{ia}])
        hold off
        legend({'Exact steady solution', 'Steady state of manufactured solution'})
        uSteady_FigureLegendStrings = [uSteady_FigureLegendStrings...
            {['$a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$']}]; %#ok<AGROW>
        %% Plot unsteady solution
        figure();
        hold on
        Times = 0:0.25:1;
        Colors = cool(length(Times));
        for it = 1:length(Times)
            plot(xh, f_um(Times(it), xh), 'Color', Colors(it,:))
        end
        hold off
        xlabel('$\tilde{x}$', 'Interpreter', 'latex');
        ylabel('$\tilde{u}$', 'Interpreter', 'latex');
        legend({'$t = 0$', '$t = 0.25$', '$t = 0.5$', '$t = 0.75$', '$t = 1$'},...
            'Interpreter', 'latex')
        title(['$a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$'])
        %% Plot derived Neumann boundary condition
        figure(h_Figure)
        hold on
        th = 0:0.01:1;
        plot(th, f_hm(th), [PeStyles{iPe}, aColors{ia}])
        hold off
        h_FigureLegendStrings = [h_FigureLegendStrings,...
            {['$a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$']}]; %#ok<AGROW>
        %% Plot derived source term
        figure()
        hold on
        for it = 1:length(Times)
            plot(xh, f_sm(Times(it), xh), 'Color', Colors(it,:))
        end
        hold off
        xlabel('$\tilde{x}$', 'Interpreter', 'latex');
        ylabel('$\tilde{s}$', 'Interpreter', 'latex');
        legend({'$t = 0$', '$t = 0.25$', '$t = 0.5$', '$t = 0.75$', '$t = 1$'},...
            'Interpreter', 'latex')
        title(['$a = ', sprintf('%g', a), ...
            '$, $\mathrm{Pe}_r = ', sprintf('%g', Pe_r), '$'])
    end
end

figure(h_Figure)
legend(h_FigureLegendStrings,  'Interpreter', 'latex')

figure(uSteady_Figure)
legend(uSteady_FigureLegendStrings,  'Interpreter', 'latex')