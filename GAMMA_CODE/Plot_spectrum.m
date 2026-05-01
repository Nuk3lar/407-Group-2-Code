%% Parameters
element = "Eu-152";
materials = { "W-PLA High Density", "Cu-PLA", "Graded-1", "Graded-2", "Graded-3"};
thicknesses = { ["3", "4"], ["1", "3"], ["5"] , ["5"] , ["5"] };

%% --- Colorblind-friendly palette (Okabe-Ito) ---
cbColors = [...
    86/255 180/255 233/255;  % sky blue
    213/255 94/255 0;        % vermillion
    230/255 159/255 0;       % orange
    0 158/255 115/255;       % bluish green
    0 114/255 178/255;       % blue
    204/255 121/255 167/255; % reddish purple
    
    % --- additional safe extensions ---
    0.35 0.35 0.35;          % dark grey
    0.6 0.6 0.6;             % light grey
    0.8 0.47 0.65;           % muted pink
    0.4 0.65 0.85];          % soft blue

colorIdx = 1; % index for cycling through colors

%% Initialize figure
figure('Units','centimeters', 'Position',[5 5 12 12], Theme="light");
hold on

plotHandles = gobjects(0);
legendLabels = strings(0);

%% --- Define graded shields ---
gradedLibrary = containers.Map;

gradedLibrary("Graded-1") = struct( ...
    'materials', ["W-PLA High Density", "Cu-PLA"], ...
    'thicknesses_mm', ["3","1"]);

gradedLibrary("Graded-2") = struct( ...
    'materials', ["W-PLA High Density","Cu-PLA"], ...
    'thicknesses_mm', ["3","1.5"]);

gradedLibrary("Graded-3") = struct( ...
    'materials', ["W-PLA High Density","Cu-PLA"], ...
    'thicknesses_mm', ["3.5","1.5"]);

%% --- Plot baseline (BLACK) ---
baselineFile = sprintf('Spectrum Data//%s Baseline.csv', element);
Tbase = readtable(baselineFile);
xBase = Tbase.Energy;
yBase = Tbase.Data;
yBase(yBase < 0) = 0;

hBase = plot(xBase, yBase, 'LineWidth', 1.5, 'Color', [0 0 0]); % BLACK
plotHandles(end+1) = hBase;
legendLabels(end+1) = sprintf('%s Baseline', element);

%% --- Loop over materials ---
for m = 1:length(materials)
    matName = materials{m};
    matThicknesses = thicknesses{m};
    
    % --- CASE 1: Graded material ---
    if startsWith(matName, "Graded")
        
        shieldDef = gradedLibrary(matName);
        
        % Build filename string first
        fileStr = sprintf('%s', element);
        
        for g = 1:length(shieldDef.materials)
            mat_g   = shieldDef.materials(g);
            thick_g = shieldDef.thicknesses_mm(g);
            
            fileStr = fileStr + ...
                        sprintf(" %smm %s", thick_g, mat_g);
        end
        
        % Use fullfile here
        filename = fullfile(pwd, 'Spectrum Data', fileStr + ".csv");
        
        % Read and plot
        T = readtable(filename);
        x = T.Energy;
        y = T.Data;
        y(y < 0) = 0;
        
        currentColor = cbColors(mod(colorIdx-1, size(cbColors,1)) + 1, :);
        colorIdx = colorIdx + 1;
        
        h = plot(x, y, 'LineWidth', 1.5, 'Color', currentColor);
        plotHandles(end+1) = h;
        
        legendLabels(end+1) = matName;
        
    else
        % --- CASE 2: Normal material ---
        for i = 1:length(matThicknesses)
            
            fileStr = sprintf('%s %smm %s.csv', ...
                element, matThicknesses(i), matName);
            
            % Use fullfile here
            filename = fullfile(pwd, 'Spectrum Data', fileStr);
            
            T = readtable(filename);
            x = T.Energy;
            y = T.Data;
            y(y < 0) = 0;
            
            currentColor = cbColors(mod(colorIdx-1, size(cbColors,1)) + 1, :);
            colorIdx = colorIdx + 1;
            
            h = plot(x, y, 'LineWidth', 1.5, 'Color', currentColor);
            plotHandles(end+1) = h;
            
            legendLabels(end+1) = sprintf('%s %smm', matName, matThicknesses(i));
        end
    end
end

%% Labels and formatting
xlabel('Energy (keV)');
ylabel('Counts Per Channel');
set(gca, 'YScale', 'log')
xlim([0 400])
ylim([10^(3) 1e5])

legend(plotHandles, legendLabels, 'Location','northeast');
grid on

%% Save figure
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12 12]);
filename = sprintf('lab_spectrum_%s.png', element);
exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'none');