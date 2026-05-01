%% Parameters
energy = 180;       % kV for the CSV files
nPhsp = 3;      % number of filters
span = 0.05;       % Loess smoothing fraction (not used with Savitzky-Golay here)
filter = 1;     % Filter Used in plotting

%% Initialize figure (10x10 cm)
figure('Units','centimeters', 'Position',[5 5 12 12]); % [left bottom width height]
hold on

plotHandles = gobjects(nPhsp,1);   % store plot handles for legend
legendLabels = strings(nPhsp,1);   % store legend labels
x_scaled = [];                         % placeholder for X-axis

%% Loop through CSV files
for i = 1:nPhsp
    % Build file name dynamically
    filename = sprintf('XStrahl Data//%dkvP_Filter_%d_%d.csv', energy, filter, i);
    
    % Read CSV
    T = readtable(filename);
    
    % Scale X-axis (assumes all CSVs have same Var1)
    if isempty(x_scaled)
        x_scaled = T.Var1 * 1000;
    end
    
    % Get Y data
    y = T.Var2;
    
    % Apply Savitzky-Golay smoothing
    windowSize = 101;  % adjust for peak width
    polyOrder = 4;
    y_smooth = sgolayfilt(y, polyOrder, windowSize);
    
    % Clip negatives
    y_smooth(y_smooth < 0) = 0;
    
    % Plot smoothed data
    plotHandles(i) = plot(x_scaled, y_smooth, 'LineWidth',1.5);
    
    % Prepare legend label
    switch(i)
        case 1
            legendLabels(i) = sprintf('Before Filter');

        case 2
            legendLabels(i) = sprintf('After Filter');

        case 3
            legendLabels(i) = sprintf('End of Applicator');

    end
    
end

%% Labels, title, and legend
xlabel('Energy (keV)');
ylabel('$\mathrm{Planar\ fluence}\;\mathrm{[Counts/(incident\ particle\cdot cm^{2}\cdot keV)]}$', ...
       'Interpreter','latex');
title(sprintf('%dkV Planar Fluence', energy));
legend(plotHandles, legendLabels, 'Location','best');
grid on
linkdata on;

%% Optional: save figure at 300 DPI for publication
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 12]);  % 12x12 cm
% Construct filename and save at 300 DPI
filename = sprintf('beam_spectrum_%dkvP_filter_areas.png', energy);
print(gcf, filename, '-dpng', '-r300');
