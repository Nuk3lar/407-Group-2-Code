%% Parameters
energies = [80 120 160 200 300];   % energies to plot (kV)
filter   = 1;              % always Filter_1
phsp     = 2;              % ONLY plot PHSP2

%% Initialize figure (12x12 cm)
figure('Units','centimeters', 'Position',[5 5 24 12]);
hold on

plotHandles  = gobjects(numel(energies),1);
legendLabels = strings(numel(energies),1);

xMaxAll = 0;   % track global X max

%% Loop over energies
for e = 1:numel(energies)

    energy = energies(e);

    % Read CSV
    filename = sprintf('XStrahl Data/%dkvP_Filter_%d_%d.csv', ...
                        energy, filter, phsp);
    T = readtable(filename);

    % X and Y (each dataset has its own X!)
    x = T.Var1 * 1000;   % keV
    y = T.Var2;

    % Track global X max
    xMaxAll = max(xMaxAll, max(x));

    % Savitzky–Golay smoothing
    windowSize = 101;
    polyOrder  = 4;
    y_smooth = sgolayfilt(y, polyOrder, windowSize);

    % Clip negatives
    y_smooth(y_smooth < 0) = 0;

    % Plot using *its own X*
    plotHandles(e) = plot(x, y_smooth, 'LineWidth', 1.5);

    legendLabels(e) = sprintf('%dkV', energy);
end

%% Axes, labels, legend
xlabel('Energy (keV)');
ylabel('$\mathrm{Planar\ fluence}\;\mathrm{[Counts/(incident\ particle\cdot cm^{2}\cdot keV)]}$', ...
       'Interpreter','latex');
title('Planar Fluence at PHSP2 (Filter 1)');
legend(plotHandles, legendLabels, 'Location','best');
grid on

% Enforce common axes
xlim([0 xMaxAll]);

linkdata on

%% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 24 12]);
print(gcf, 'beam_spectrum_PHSP2_all_energies.png', '-dpng', '-r300');
