syms x Per a g h u(x)
ue = dsolve([(a*diff(u, x) - 1/Per*diff(u, x, 2) == 0),...
    (-1/Per*subs(diff(u, x), x, 0) == h),...
    (subs(u, x, 1) == g)]);
display(ue)
%%
ue_a0 = limit(ue, a, 0);
display(ue_a0)
%%
h0 = solve(subs(ue, x, 0));
display(h0)
%%
h0_a0 = limit(h0, a, 0);
display(h0_a0)
%%
ue_h0 =  subs(ue, h, h0);
ue_h0 = simplify(ue_h0);
display(ue_h0);
%%
ue_a0_h0 = subs(ue_a0, h, h0_a0);
ue_a0_h0 = simplify(ue_a0_h0);
display(ue_a0_h0);
%%
xh = 0:0.01:1;
gList = [-1, -2];
PerList = [1, 0.5];
Styles = {'-', '--'};
aList = [0, -1, -5];
Colors = {'k', 'b', 'r'};
for ig = 1:length(gList)
    g = gList(ig);
    figure()
    title(['g = ', sprintf('%g', g)])
    hold on
    LegendStrings = {};
    for iPer = 1:length(PerList)
        Per = PerList(iPer);
        Style = Styles{iPer};
        for ia = 1:length(aList)
            a = aList(ia);
            Color = Colors{ia};
            if a == 0
                f = matlabFunction(ue_a0_h0);
                uh = f(g, xh);
            else
                f = matlabFunction(ue_h0);
                uh = f(Per, a, g, xh);
            end
            plot(xh, uh, [Style, Color])
            LegendStrings = [LegendStrings,...
                {['Pe_r = ', sprintf('%g', Per),...
                ', a = ', sprintf('%g', a)]}]; %#ok<AGROW>
        end
    end
    legend(LegendStrings)
    hold off
end
