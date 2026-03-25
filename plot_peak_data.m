cd(fileparts(mfilename('fullpath')))
analyse = 3; % <-- General Purpose switch
analse_every = false;
egs_data = true;
plot_graphs = [false, false, false, true];
close all;
switch analyse
    case 1
        elements = ["Eu-152", "Cs-137"];
        materials = { {"LM-95", 'Pb'}, ...
                      {"LM-95", 'Pb'}};
        % IN MM
        thicknesses = { ...
            { {"4.21"}, {"2.68"} }, ...
            { {"4.21"}, {"2.68"} } ...
        };
        lowx = 100;
        highx = 1000;
    case 2
        elements = ["Eu-152", "Cs-137", "Am-241"];
        materials = { {"W-PLA High Density", "W-PLA Low Density", "Cu-PLA" }, ...
                      {"W-PLA High Density", "Cu-PLA"}, ...
                      {"Cu-PLA"} };
        
        % IN MM
        thicknesses = { ...
            { {"1","2","3","4"}, {"3"}, {"1.5", "3","4.5"} }, ...
            { {"1","2","3","4"}, {"1","3"}, {"1", "3"} }, ...
            { {"1", "3", "4"} } ...
        };
        lowx = 25;
        highx = 1000;
    case 3
        elements = [ "Eu-152", "Cs-137", "Am-241" ];
        materials = { {"W-PLA High Density"},...
                      {"W-PLA High Density"},...
                      {"Cu-PLA"} };
        thicknesses = { ...
            { {"1", "1.5", "2", "3", "4", "7.5"}},...
            { {"1", "2", "3", "4"} },...
            { {"1", "3", "4"} }...
        };
        lowx = 25;
        highx = 1000;    
    case 4
        elements = ["Eu-152", "Am-241"];
        materials = { {"Cu-PLA", "W-PLA High Density", "W-PLA Low Density"}, {"Cu-PLA"} };
        
        % IN MM
        thicknesses = { ...
            { {"1", "1.5", "3","4.5"}, {"1", "1.5", "2", "3", "4", "7.5"}, {"3"} }, ...
            { {"1", "3", "4"} }, ...
        };
        lowx = 25;
        highx = 1000;  
    case 5
        elements = [ "Eu-152", "Cs-137" ];
        materials = { {"W-PLA High Density", "Graded-1", "Graded-2", "Graded-3"}...
                      {"W-PLA High Density", "Graded-1"} };
        thicknesses = { ...
            { {"1", "1.5", "2", "3", "4"}, {""}, {""}, {""} },...
            { {"1", "2", "3", "4"}, {""}  }...
        };

        lowx = 100;
        highx = 1000;  
end
% --- File to store cached results ---
cacheFile = 'peak_cache.mat';

% Check if cache exists
if analse_every
    useCache = false;

elseif exist(cacheFile,'file')
    S = load(cacheFile, 'lastAnalyse', 'results', 'B_table');
    useCache = false;
    if isfield(S,'lastAnalyse') && S.lastAnalyse == analyse
        % Analyse hasn't changed use cached results
        results = S.results;
        B_table = S.B_table;
        useCache = true;
        fprintf('Using cached results for analyse = %d\n', analyse);
    end
else
    useCache = false;
end

% If no cache or analyse changed compute
if ~useCache
    fprintf('Running peak_process for analyse = %d ...\n', analyse);
    [results, B_table] = peak_process_multi(elements, materials, thicknesses);
    if egs_data
        results = [results; peak_process_egs()];
    end
    % Save cache
    lastAnalyse = analyse; 
    save(cacheFile, 'lastAnalyse', 'results', 'B_table');
end

% Graph options
graphx_cm = 12;
graphy_cm = 12;
forceZeroIntercept = true;

% --- Colorblind-friendly palette (Okabe-Ito) ---
cbColors = [...

    86/255 180/255 233/255;    % sky blue
    213/255 94/255 0;          % vermillion
    230/255 159/255 0;         % orange
    0 158/255 115/255;         % bluish green
    0 114/255 178/255;         % blue
    204/255 121/255 167/255];  % reddish purple


%PLOTTING :)
if plot_graphs(1)
    %% --- Prepare data ---
    % Remove BASELINE rows
    dataToPlot = results;
    isBaseline = dataToPlot.Thickness == "BASELINE";
    dataToPlot = dataToPlot(~isBaseline, :);
    
    % Convert Thickness from string → double
    dataToPlot.Thickness = str2double(dataToPlot.Thickness);
    
    % Create series label (Material + Thickness)
    seriesLabels = dataToPlot.Material + " " + string(dataToPlot.Thickness) + "mm";
    uniqueSeries = unique(seriesLabels);
    
    % Extract base materials
    materialsUnique = unique(dataToPlot.Material);
    baseMaterials = materialsUnique;
    baseMaterials = unique(baseMaterials);
    
    %% --- Style definitions ---
    markers = {'x', 'o', '*', 's', '^', 'v', '<', '>', 'p', 'h'};
    numThickness = length(uniqueSeries);
    numMaterials  = length(baseMaterials);
    thicknessColors = cbColors(mod(0:numThickness-1, size(cbColors,1))+1, :);
    materialColors  = cbColors(mod(0:numMaterials-1, size(cbColors,1))+1, :);
    
    figure('Units','centimeters', 'Position',[5 5 graphx_cm graphy_cm], Theme="light");
    
    hold on; box on;
    xLimits = [lowx highx];
    
    %% ==========================================================
    %% Plot ln(A/A0) vs Thickness for each Photopeak
    %% ==========================================================
    
    validPeaks = dataToPlot.PhotoPeak(~isnan(dataToPlot.PhotoPeak));
    uniquePeaks = unique(validPeaks);
    uniqueMaterials = unique(dataToPlot.Material);
    
    hold on
    
    for pIdx = 1:length(uniquePeaks)
    
        peakEnergy = uniquePeaks(pIdx);
        if isnan(peakEnergy)
            continue
        end
        for m = 1:length(uniqueMaterials)
            
            matName = uniqueMaterials(m);
    
            % Select this photopeak + material
            idx = dataToPlot.PhotoPeak == peakEnergy & ...
                  dataToPlot.Material == matName;
    
           
    
            t = dataToPlot.Thickness(idx);
            y = dataToPlot.ln_transmisson(idx);
            yErr = dataToPlot.ln_transmisson_unc(idx);
    
            % Sort by thickness (important for nice lines)
            [t, sortIdx] = sort(t);
            y = y(sortIdx);
            y(y > 0) = 0;
            yErr = yErr(sortIdx);
    
            % ---- Marker style by Photopeak ----
            markerStyle = markers{mod(pIdx-1,length(markers))+1};
    
            % ---- Line color by Material ----
            matIdx = find(baseMaterials == matName);
            colorStyle = materialColors(matIdx,:);
    
            % Plot data
            errorbar(t, y, yErr, markerStyle, ...
                'MarkerFaceColor', colorStyle, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 7, ...
                'Color', colorStyle, ...
                'LineStyle','none', ...
                'HandleVisibility','off', ...
                'DisplayName', sprintf('%s - %.0f keV', matName, peakEnergy));
            % ---- Linear Beer-Lambert Fit ----
            
            tLine = linspace(0, max(dataToPlot.Thickness)+4, 200);
            
            if sum(idx) == 1
            
                % Single point case
                if forceZeroIntercept
                    slope = y / t;        % line through origin
                    yLine = slope * tLine;
                else
                    slope = y / t;        % fallback anyway (cannot polyfit 1 point)
                    yLine = slope * tLine;
                end
            else
            
                if forceZeroIntercept
                    % ===== Constrained Fit (Intercept = 0) =====
                    % Solve least squares for y = m t
                    slope = sum(t .* y) / sum(t.^2);
                    yLine = slope * tLine;
            
                else
                    
                    % ===== Standard Unconstrained Fit =====
                    pFit = polyfit(t, y, 1);
                
                    slope = pFit(1);        % = -μ
                    intercept = pFit(2);    % = ln(B)
                
                    mu = -slope;
                    B  = exp(intercept);
                
                    yLine = polyval(pFit, tLine);
                
                    % ---- Store results ----
                    results(m, pIdx).Energy = peakEnergy;
                    results(m, pIdx).Material = matName;
                    results(m, pIdx).mu = mu;
                    results(m, pIdx).B  = B;
                   
                end
            
            end
    
    % ---- Plot fit line ----
    plot(tLine, yLine, '-', ...
        'Color', colorStyle, ...
        'LineWidth', 1, ...
        'HandleVisibility','off');
    
            % Plot fit line
            plot(tLine, yLine, '-', ...
                'Color', colorStyle, ...
                'LineWidth', 1, ...
                'HandleVisibility','off');
            
            
    
        end
        
    end
    % ============================================
    % BUILD ONE COMBINED LEGEND
    % ============================================
    
    hLegend = gobjects(0);
    labelsLegend = {};
    
    %% ----- Material color entries -----
    for m = 1:length(uniqueMaterials)
    
        matName = uniqueMaterials(m);
        matIdx = find(baseMaterials == ...
            matName);
    
        if isempty(matIdx)
            continue
        end
    
        colorStyle = materialColors(matIdx,:);
    
        h = plot(NaN, NaN, '-', ...
            'Color', colorStyle, ...
            'LineWidth', 2);
    
        hLegend(end+1) = h;
        labelsLegend{end+1} = char(matName);
    end
    
    %% ----- Spacer (optional visual separation) -----
    hSpacer = plot(NaN,NaN,'w');
    hLegend(end+1) = hSpacer;
    labelsLegend{end+1} = ' ';   % blank line
    
    %% ----- Photopeak marker entries -----
    for pIdx = 1:length(uniquePeaks)
    
        peakEnergy = uniquePeaks(pIdx);
    
        if isnan(peakEnergy)
            continue
        end
    
        markerStyle = markers{mod(pIdx-1,length(markers))+1};
    
        h = plot(NaN, NaN, markerStyle, ...
            'MarkerSize',8, ...
            'MarkerFaceColor','k', ...
            'Color','k', ...
            'LineStyle','none');
    
        hLegend(end+1) = h;
        labelsLegend{end+1} = sprintf('%.0f keV', peakEnergy);
    end
    
    % ---- Proper legend spacer ----
    legend_spacer = plot(NaN, NaN, '-', ...
        'Color', 'w', ...          % white (invisible on white bg)
        'LineWidth', 0.01);
    
    hLegend(end+1) = legend_spacer;
    labelsLegend{end+1} = ' ';     % blank label
    
    %% ----- HVL & TVL -----
    
    hHVL = plot(NaN,NaN,'--k','LineWidth',1.2);
    hLegend(end+1) = hHVL;
    labelsLegend{end+1} = 'HVL $(A/A_0 = 0.5)$';
    yHVL = log(0.5);
    
    yline(yHVL, '--k', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'HVL (A/A_0 = 0.5)', 'Interpreter', 'latex');

    %{
    hTVL = plot(NaN,NaN,'-k','LineWidth',1.2);
    hLegend(end+1) = hTVL;
    labelsLegend{end+1} = 'TVL $(A/A_0 = 0.1)$';

    yTVL = log(0.1);
    yline(yTVL, '-k', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'TVL (A/A_0 = 0.1)', 'Interpreter', 'latex');
    

    hCVL = plot(NaN,NaN,':k','LineWidth',1.2);
    hLegend(end+1) = hCVL;
    labelsLegend{end+1} = 'CVL $(A/A_0 = 0.01)$';
  
    yCVL = log(0.01);
    yline(yCVL, ':k', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'CVL $(A/A_0 = 0.01)$', 'Interpreter', 'latex');
    %}
    %% ----- Create final legend -----
    legend(hLegend, labelsLegend, ...
        'Location','eastoutside', ...
        'Interpreter','latex', ...
        'FontSize',9);
    
    set(gca,'XScale','linear','YScale','linear');
    xlabel('Thickness (mm)', ...
           'FontSize',12, ...
           'FontWeight','bold', ...
           'Interpreter','latex');
    xlim([0 max(dataToPlot.Thickness)+2.5]);
    ylim([ -1 0 ]);
    ylabel('$\ln(A/A_0)$', ...
           'FontSize',12, ...
           'FontWeight','bold', ...
           'Interpreter','latex');
    grid on
    
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 graphx_cm graphy_cm]);  % 12x12 cm
    filename = sprintf('thickness_attenuation_intercept_%d.png', forceZeroIntercept);
    exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'white')
end



% ----------------------------------------------------------
if plot_graphs(2)
    
    %% --- Prepare data ---
    dataToPlot = results(results.Material ~= "BASELINE", :);
    
    % Create series label (Material + Thickness)
    seriesLabels = dataToPlot.Material + " " + string(dataToPlot.Thickness) + "mm";
    uniqueSeries = unique(seriesLabels);
    uniqueThickness = unique(dataToPlot.Thickness);
    
    % Extract base materials
    materialsUnique = unique(dataToPlot.Material);
    baseMaterials = erase(materialsUnique, [" High Density"," Low Density"]);
    baseMaterials = unique(baseMaterials);
    
    %% --- Style definitions ---
    markers = {'<', '^', 'v', '>', 'o', 'p', 'h', 's'};
    

    
    materialColors  = cbColors(mod(0:numMaterials-1, size(cbColors,1))+1, :);    
    figure('Units','centimeters', 'Position',[5 5 graphx_cm graphy_cm], Theme="light");
    
    hold on; box on;
    
    %% ==========================================================
    %% 1) Plot measured data (marker by material, colour by thickness)
    %% ==========================================================
    
    for s = 1:length(uniqueSeries)
    
        label = uniqueSeries(s);
        idx = seriesLabels == label;
    
        % Extract base material
        matName = dataToPlot.Material(find(idx,1));
        matBase = erase(matName, [" High Density"," Low Density"]);
    
        % Find material colour
        matIdx = find(baseMaterials == matBase);
        colorStyle = materialColors(matIdx,:);
    
        % Find thickness marker
        thicknessVal = dataToPlot.Thickness(find(idx,1));
        tIdx = find(uniqueThickness == thicknessVal);
        markerStyle = markers{mod(tIdx-1,length(markers))+1};
    
        x = dataToPlot.Centroid(idx);
        y = dataToPlot.mu_mass(idx);
        yErr = dataToPlot.mu_mass_unc(idx);
    
        epsVal = 1e-8;
        yLower = max(y - yErr, epsVal);
        yErrLower = y - yLower;
        yErrUpper = yErr;
    
        errorbar(x, y, yErrLower, yErrUpper, markerStyle, ...
            'MarkerFaceColor', colorStyle, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 6, ...
            'Color', colorStyle, ...
            'LineStyle','none', ...
            'LineWidth',1, ...
            'HandleVisibility','off');
    end
    
    %% ==========================================================
    %% 2) Plot XCOM curves (colour by material only)
    %% ==========================================================
    
    for m = 1:length(baseMaterials)
    
        matBase = baseMaterials(m);
        % --- Skip graded materials ---
        if startsWith(matBase, "Graded-")
            continue;
        end

        if startsWith(matBase, "(EGS)")
            continue;
        end
        colorStyle = materialColors(m,:);
    
        filename = sprintf("XCOM Data//XCOM_%s.csv", matBase);
    
        opts = detectImportOptions(filename,'NumHeaderLines',6);
        xcomTable = readtable(filename, opts);
    
        energy_keV = xcomTable{:,1} * 1e3;
        mu_xcom    = xcomTable{:,3};
    
        plot(energy_keV, mu_xcom, '-', ...
             'LineWidth',1.6, ...
             'Color',colorStyle, ...
             'HandleVisibility','off');
    end
    
    %% ==========================================================
    %% 3) Log-log cubic fit + shaded region (per material)
    %% ==========================================================
    
    
    xLimits = [lowx highx];
    
    for m = 1:length(baseMaterials)
        matBase = baseMaterials(m);
        colorStyle = materialColors(m,:);
    
        % --- Select data for this material, ignore baseline ---
        idxMat = contains(dataToPlot.Material, matBase) & dataToPlot.Material ~= "BASELINE";
        x = dataToPlot.Centroid(idxMat);
        y = dataToPlot.mu_mass(idxMat);
        yErr = dataToPlot.mu_mass_unc(idxMat);
    
        % Only positive y values
        posIdx = y > 0;
        xFit = x(posIdx);
        yFit = y(posIdx);
        yErrFit = yErr(posIdx);
    
        % --- Compute mean relative uncertainty (fraction) ---
        relError = 10 * mean(yErrFit ./ yFit, 'omitnan');  % e.g., 0.15
        relErrorPct = relError * 100;                % 15 for 15%
    
        % --- Cubic fit in log-log space ---
        p = polyfit(log(xFit), log(yFit), 3);
    
        xLine = logspace(log10(xLimits(1)), log10(xLimits(2)), 4000);
        yLine = exp(polyval(p, log(xLine)));
        
        % Shaded uncertainty region
        yUpper = yLine * (1 + relError);
        yLower = yLine * (1 - relError);
    
        % Plot shaded region with legend label including %
        fill([xLine fliplr(xLine)], ...
             [yUpper fliplr(yLower)], ...
             colorStyle, ...
             'FaceAlpha',0.15, ...
             'EdgeColor','k', ...
             'HandleVisibility','off'); 
        
        % Plot fit line
        plot(xLine, yLine, '--', ...
             'LineWidth',1.6, ...
             'Color',colorStyle, ...
             'HandleVisibility','off');
             %'DisplayName', sprintf('%s Fit $\\pm$ $10 \\cdot\\overline{\\delta\\mu_m}$ (%.1f\\%%) ', matBase, (relErrorPct / 10)));
    end
    
    %% ==========================================================
    %% 4) Final axis formatting
    %% ==========================================================
    
    set(gca,'XScale','log','YScale','log');
    xlim([lowx highx]);
    ylim([1e-2 1e2]);
    xlabel('Photon Energy (keV)', ...
           'FontSize',12, ...
           'Interpreter','latex');
    ylabel('Mass Attenuation Coefficient $\mu/\rho$ (cm$^2$/g)', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    %{
    title('Mass Attenuation (\mu/\rho) with XCOM Reference', ...
           'FontSize',12);
    %}
    grid on;
    
    %% ============================================
    %% CUSTOM LEGEND FIXED
    %% ============================================
    
    hLegend = gobjects(0);
    labelsLegend = {};
    
    % --- MATERIALS (solid lines) ---
    for m = 1:length(baseMaterials)
        colorStyle = materialColors(m,:);
        % create a line handle that matches XCOM data line
        h = plot(NaN, NaN, '-', ...
                 'Color', colorStyle, ...
                 'LineWidth', 1.8);
        hLegend(end+1) = h;
        labelsLegend{end+1} = baseMaterials(m);
    end
    
    % Spacer
    hLegend(end+1) = plot(NaN, NaN, 'w');
    labelsLegend{end+1} = ' ';
    
    % Shaded uncertainty region
    hShade = patch(NaN, NaN, [0.7 0.7 0.7], ...
                   'FaceAlpha', 0.15, ...
                   'EdgeColor', 'k');
    hLegend(end+1) = hShade;
    labelsLegend{end+1} = '$\mu_m \pm 10 \times \overline{\delta\mu_m}$';
    
    % XCOM data (dashed line)
    hXCOM = plot(NaN, NaN, '--k', 'LineWidth', 1.6);
    hLegend(end+1) = hXCOM;
    labelsLegend{end+1} = 'Fit Curve';
    
    % --- THICKNESS markers ---
    for tIdx = 1:length(uniqueThickness)
        markerStyle = markers{mod(tIdx-1, length(markers))+1};
        % use black marker with face color to show thickness
        hMarker = plot(NaN, NaN, markerStyle, ...
                       'MarkerFaceColor', 'k', ...
                       'MarkerEdgeColor', 'k', ...
                       'MarkerSize', 7, ...
                       'LineStyle', 'none');
        hLegend(end+1) = hMarker;
        labelsLegend{end+1} = sprintf('%.2f mm', uniqueThickness(tIdx));
    end
    
    % Final legend
    legend(hLegend, labelsLegend, ...
           'Location', 'eastoutside', ...
           'Interpreter', 'latex', ...
           'FontSize', 9, ...
           'Box', 'on');
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 graphx_cm graphy_cm]);  % 12x12 cm
    filename = sprintf('mass_attenuation.png');
    exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'white')
    
end


if plot_graphs(3)
    
    %% ===============================
    %% Plot B(E,x) vs Energy per material & thickness
    %% ===============================
    
    % --- Define styles ---
    markers = {'o', 's', '^', '<',  'd', 'h', 'p'};       % markers for thickness
    lineStyles = {'-', '--', ':', '-.'};     % lines for thickness
    % --- Assign color to each material ---
    materialsFull = unique(B_table.Material);
    numMaterials  = length(materialsFull);
    
    % Cycle through cbColors if you have more materials than colors
    materialColors = cbColors(mod(0:numMaterials-1, size(cbColors,1)) + 1, :);
    
    % Convert Thickness in results table to numeric, BASELINE → NaN
    thickVals = nan(height(results),1);
    for i = 1:height(results)
        val = results.Thickness{i};
        if ~strcmp(val,"BASELINE")
            thickVals(i) = str2double(val);
        end
    end
    uniqueThickness = unique(thickVals(~isnan(thickVals)));
    
    figure('Units','centimeters', 'Position',[5 5 graphx_cm graphy_cm], Theme="light");
    hold on; box on;
    xlabel('Energy (keV)', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    ylabel('$B_{(E,x)}$', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    
    %% --- Loop over materials and thicknesses ---
    for m = 1:length(materialsFull)
        matName = materialsFull{m};
        colorStyle = materialColors(m,:);
    
        % --- Optional: Plot B(E,0) from BASELINE ---
        idxBaseline = strcmp(B_table.Material,"BASELINE");
        if any(idxBaseline)
            E_baseline = B_table.Energy(idxBaseline);
            B0_mean   = B_table.B0_mean(idxBaseline);
            B0_delta  = B_table.B0_delta(idxBaseline);
    
            errorbar(E_baseline, B0_mean, B0_delta, 's', ...
                'MarkerFaceColor', [0 0 0], ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 6, ...
                'Color', [0 0 0], ...
                'LineStyle','--', ...
                'LineWidth', 1.5, ...
                'DisplayName', 'BASELINE B(E,0)');
        end
    
        for t = 1:length(uniqueThickness)
            thick = uniqueThickness(t);
    
            % --- Find rows matching this material & approximate thickness ---
            idxB = strcmp(B_table.Material, matName) & abs(B_table.Ideal_Thickness - thick) < 1e-6;
            if sum(idxB) == 0
                continue
            end
    
            E = B_table.Energy(idxB);
            Bmean = B_table.B_mean(idxB);
            Bdelta = B_table.B_delta(idxB);
    
            % --- Marker and line style ---
            markerStyle = markers{mod(t-1,length(markers))+1};
            lineStyle   = lineStyles{mod(t-1,length(lineStyles))+1};
    
            % --- Plot data points ---
            errorbar(E, Bmean, Bdelta, markerStyle, ...
                'MarkerFaceColor', colorStyle, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 6, ...
                'Color', colorStyle, ...
                'LineStyle','none', ...
                'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s, %.1f mm', matName, thick));
    
            % --- Fit polynomial with forced (0,1) point ---
            nPts = length(E);
            deg = min(3, nPts);  % include forced point
            if deg >= 1
                % Append forced point at E=0, B=1
                EfitData = [0; E];
                BfitData = [1; Bmean];
    
                % Fit polynomial
                pFit = polyfit(EfitData, BfitData, deg);
    
                % Evaluate fit for plotting
                Eplot = linspace(0, highx, 100);
                Bplot = polyval(pFit, Eplot);
    
                % Plot polynomial fit
                plot(Eplot, Bplot, lineStyle, 'Color', colorStyle, 'LineWidth', 1.5, 'HandleVisibility','off');
            end
        end
    end

    %% ============================================
    %% CUSTOM LEGEND
    %% ============================================
    
    hLegend = gobjects(0);
    labelsLegend = {};
    
    % --- 1) Material colors ---
    materialsFull = unique(B_table.Material);
    materialColors = cbColors(mod(0:numMaterials-1, size(cbColors,1)) + 1, :);
    
    for m = 1:length(materialsFull)
        matName = materialsFull{m};
        colorStyle = materialColors(m,:);
        
        % Solid line as placeholder for material
        h = plot(NaN, NaN, '-', 'Color', colorStyle, 'LineWidth', 1.8);
        hLegend(end+1) = h;
        labelsLegend{end+1} = sprintf('%s', matName);
    end
    
    % Spacer
    spacer = plot(NaN, NaN, 'w');
    hLegend(end+1) = spacer;
    labelsLegend{end+1} = ' ';
    
    % --- 2) Thickness markers ---
    uniqueThickness = unique(B_table.Ideal_Thickness(~isnan(B_table.Ideal_Thickness)));
    
    for t = 1:length(uniqueThickness)
        thick = uniqueThickness(t);
        
        markerStyle = markers{mod(t-1,length(markers))+1};
        lineStyle   = lineStyles{mod(t-1,length(lineStyles))+1};
        
        % Plot dummy marker/line for legend
        h = plot(NaN, NaN, lineStyle, 'Color', [0 0 0], 'Marker', markerStyle, ...
            'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hLegend(end+1) = h;
        labelsLegend{end+1} = sprintf('%.1f mm', thick);
    end
    % Spacer
    spacer = plot(NaN, NaN, 'w');
    hLegend(end+1) = spacer;
    labelsLegend{end+1} = ' ';
    ybe0 = 1;
    
    yline(ybe0, '-k', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'HVL (A/A_0 = 0.5)', 'Interpreter', 'latex');
    be0 = plot(NaN,NaN,'-k','LineWidth',1.2);
    hLegend(end+1) = be0;
    labelsLegend{end+1} = '$B_{(E, 0)}$ (1)';

    % --- Create legend ---
    legend(hLegend, labelsLegend, ...
           'Location','eastoutside', ...
           'Interpreter','latex', ...
           'FontSize',9, ...
           'Box','on');
   set(gca, 'YScale', 'linear');   % linear Y axis
    ylim([0.8 1.6]);                    % fix limits
    xlim([0 highx]);
    grid on;
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 graphx_cm graphy_cm]);  % 12x12 cm
    filename = sprintf('broad_beam.png');
    exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'white')
end
warning('off','MATLAB:xlswrite:AddSheet'); %optional
writetable(results,'results.xlsx','Sheet',1);





if plot_graphs(4)


    T = results;
    tol = 2;
    T(T.mu_linear <= 0, :) = [];
   

    % Remove baseline + NaNs
    T(strcmp(T.Material, "BASELINE"), :) = [];
    T(isnan(T.mu_linear), :) = [];

    % Merge materials
    T.Material(contains(T.Material, "W-PLA")) = "W-PLA";
    T.Material(contains(T.Material, "Cu-PLA")) = "Cu-PLA";

    % Convert thickness
    T.thick_num = str2double(T.Thickness);

    % Split datasets
    T_exp = T(T.Monte_Carlo == false, :);
    T_mc  = T(T.Monte_Carlo == true,  :);

    %% --- Match EXP ↔ MC using tolerance ---
    mu_exp = [];
    mu_exp_delta = [];
    mu_mc  = [];
    mu_mc_delta  = [];
    thick  = [];
    mat    = strings(0,1);
    energy = [];

    for i = 1:height(T_exp)

        row = T_exp(i,:);
        fprintf("Element: %s, Mat: %s, Peak %dkeV\n", row.Element, row.Material, row.PhotoPeak)
        % Match element/material/thickness (numeric!)
        idx = strcmp(T_mc.Element, row.Element) & ...
              strcmp(T_mc.Material, row.Material) & ...
              abs(T_mc.thick_num - row.thick_num) < 1e-6;

        candidates = T_mc(idx,:);

        if isempty(candidates)
            continue
        end

        % Closest energy
        [dmin, ind] = min(abs(candidates.PhotoPeak - row.PhotoPeak));

        if dmin < tol
            mu_exp(end+1,1) = row.mu_linear;
            mu_exp_delta(end+1,1) = row.mu_unc;
            mu_mc(end+1,1)  = candidates.mu_linear(ind);
            mu_mc_delta(end+1,1)  = candidates.mu_unc(ind);
            thick(end+1,1)  = row.thick_num;
            mat(end+1,1)    = row.Material;
            energy(end+1,1) = row.PhotoPeak;
        end
    end

    %% --- Compute ratio ---
    ratio = mu_exp ./ mu_mc;
    ratio_unc = ratio .* sqrt( ...
        (mu_exp_delta ./ mu_exp).^2 + ...
        (mu_mc_delta ./ mu_mc).^2 ...
    );
    
    % Clean invalid values
    ratio_unc(~isfinite(ratio_unc)) = NaN;

    if isempty(ratio)
        error('No matched data found — check tolerance or inputs.');
    end

    %% --- Auto-cluster energies ---
    energy_sorted = sort(energy);
    groups = {};
    current_group = energy_sorted(1);

    for i = 2:length(energy_sorted)
        if abs(energy_sorted(i) - current_group(end)) < tol
            current_group(end+1) = energy_sorted(i);
        else
            groups{end+1} = current_group;
            current_group = energy_sorted(i);
        end
    end
    groups{end+1} = current_group;

    % Representative energy per group

    %% --- Combined Legend for Materials + Photopeaks ---
    figure('Units','centimeters', 'Position',[5 5 graphx_cm graphy_cm], Theme="light");
    hold on;
    
    % --- Define markers and line styles ---
    markers    = {'o', 's', '^', 'h', 'p', '*', 'v', '>'};
    lineStyles = {'-', '--', ':', '-.'};
    
    materials      = unique(mat);
    energy_centers = unique(energy);
    colors         = cbColors(mod(0:length(materials)-1, size(cbColors,1))+1, :);
    
    % --- Plot actual data ---
    for i = 1:numel(materials)
        for j = 1:numel(energy_centers)
    
            idx = strcmp(mat, materials(i)) & abs(energy - energy_centers(j)) < tol;
    
            if any(idx)
                % Sort by thickness
                [th_sorted, order] = sort(thick(idx));
                r_sorted    = ratio(idx);
                r_sorted    = r_sorted(order);
                r_unc_sorted = ratio_unc(idx);
                r_unc_sorted = r_unc_sorted(order);
    
                % --- Style ---
                lineStyle = lineStyles{mod(j-1,numel(lineStyles))+1};
                marker    = markers{mod(j-1,numel(markers))+1};
    
                % --- Errorbar plot ---
                errorbar(th_sorted, r_sorted, r_unc_sorted, ...
                    'Color', colors(i,:), ...
                    'LineStyle', lineStyle, ...
                    'LineWidth', 1.2, ...
                    'Marker', marker, ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', colors(i,:), ...
                    'MarkerEdgeColor', colors(i,:));
            end
        end
    end
    
    % --- Build legend: first materials (lines) ---
    hLegend = gobjects(0);
    labelsLegend = {};
    
    for i = 1:numel(materials)
        h = plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        hLegend(end+1) = h;
        labelsLegend{end+1} = char(materials(i));
    end
    
    % --- Add spacer ---
    hSpacer = plot(NaN, NaN, 'w');
    hLegend(end+1) = hSpacer;
    labelsLegend{end+1} = ' ';
    
    % --- Add photopeak markers ---
    for j = 1:numel(energy_centers)
        marker = markers{mod(j-1,numel(markers))+1};
        h = plot(NaN, NaN, marker, 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
        hLegend(end+1) = h;
        labelsLegend{end+1} = sprintf('%.0f keV', energy_centers(j));
    end
    
    legend(hLegend, labelsLegend, ...
               'Location', 'eastoutside', ...
               'Interpreter', 'latex', ...
               'FontSize', 9, ...
               'Box', 'on');
    ylabel('Experimental/EGS Linear Attenuation Ratio, $\mu_\mathrm{Exp} / \mu_\mathrm{MC}$', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    xlabel('Thickness (mm)', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    grid on;
    box on;
    xlim([0 8]);

    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 graphx_cm graphy_cm]);  % 12x12 cm
    filename = sprintf('mu_linear_ratio.png');
    exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'white')
    %% --- Prepare CSV output ---

    % Create table
    T_out = table();
    T_out.Material   = mat;
    T_out.Thickness  = thick;
    T_out.PhotoPeak  = energy;
    T_out.mu_exp     = mu_exp;
    T_out.mu_mc      = mu_mc;
    T_out.mu_ratio   = ratio;
    
    % Optional: sort by material, then peak energy, then thickness
    T_out = sortrows(T_out, {'Material','PhotoPeak','Thickness'});
    
    % Write CSV
    writetable(T_out, 'mu_ratios_by_material.csv');
    
    fprintf('CSV file written: mu_ratios_by_material.csv\n');


    %% --- μ_exp vs μ_MC plot with improved legend ---
    
    figure('Units','centimeters', 'Position',[5 5 graphx_cm graphy_cm], Theme="light");
    hold on;
    
    materials = unique(mat);
    energy_peaks = unique(energy);
    colors = cbColors(mod(0:length(materials)-1, size(cbColors,1))+1, :);
    
    markers    = {'o', 's', '^', 'h', 'p', '*', 'v', '>'};
    
    % --- Plot data ---
    for i = 1:numel(materials)
        for j = 1:numel(energy_peaks)
    
            idx = strcmp(mat, materials(i)) & abs(energy - energy_peaks(j)) < tol;
    
            if any(idx)
                [th_sorted, order] = sort(thick(idx));
    
                mu_exp_sorted = mu_exp(idx); 
                mu_exp_sorted = mu_exp_sorted(order);
                mu_exp_unc    = mu_exp_delta(idx);
                mu_exp_unc    = mu_exp_unc(order);
    
                mu_mc_sorted  = mu_mc(idx);  
                mu_mc_sorted  = mu_mc_sorted(order);
                mu_mc_unc     = mu_mc_delta(idx);
                mu_mc_unc     = mu_mc_unc(order);
    
                marker = markers{mod(j-1,numel(markers))+1};
    
                % ===== Experimental (solid, filled) =====
                errorbar(th_sorted, mu_exp_sorted, mu_exp_unc, ...
                    'Color', colors(i,:), ...
                    'LineStyle','-', ...
                    'LineWidth',1.5, ...
                    'Marker', marker, ...
                    'MarkerSize',6, ...
                    'MarkerFaceColor', colors(i,:), ...
                    'MarkerEdgeColor', colors(i,:));
    
                % ===== Simulation (dashed, open) =====
                errorbar(th_sorted, mu_mc_sorted, mu_mc_unc, ...
                    'Color', colors(i,:), ...
                    'LineStyle','--', ...
                    'LineWidth',1.5, ...
                    'Marker', marker, ...
                    'MarkerSize',6, ...
                    'MarkerFaceColor','none', ...
                    'MarkerEdgeColor', colors(i,:));
            end
        end
    end
    
    ylabel('Linear Attenuation Coefficient, $\mu$ ($\mathrm{cm}^{-1}$)', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    xlabel('Thickness (mm)', ...
           'FontSize', 12, ...
           'Interpreter', 'latex', ...
           'FontWeight','bold');
    grid on; box on;
    
    %% --- Build enhanced legend ---
    hLegend = gobjects(0);
    labelsLegend = {};
    
    % 1️⃣ Line styles explanation
    h_exp = plot(NaN, NaN, 'k-', 'LineWidth', 2.5);   % solid line = μ_exp
    h_mc  = plot(NaN, NaN, 'k--','LineWidth', 2.5);   % dashed line = μ_MC
    hLegend = [hLegend, h_exp, h_mc];
    labelsLegend = [labelsLegend, {'$\mu_\mathrm{exp}$', '$\mu_\mathrm{MC}$'}];

    % 2️⃣ Spacer
    hSpacer = plot(NaN, NaN, 'w');
    hLegend(end+1) = hSpacer;
    labelsLegend{end+1} = ' ';
    
    % 3️⃣ Materials (color lines)
    for i = 1:numel(materials)
        h = plot(NaN, NaN, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
        hLegend(end+1) = h;
        labelsLegend{end+1} = char(materials(i));
    end
    
    % 4️⃣ Spacer
    hSpacer2 = plot(NaN, NaN, 'w');
    hLegend(end+1) = hSpacer2;
    labelsLegend{end+1} = ' ';
    
    % 5️⃣ Photopeak markers
    for j = 1:numel(energy_peaks)
        marker = markers{mod(j-1,numel(markers))+1};
        h = plot(NaN, NaN, marker, 'MarkerSize',8, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
        hLegend(end+1) = h;
        labelsLegend{end+1} = sprintf('%.0f keV', energy_peaks(j));
    end
    set(gca,'XScale','linear','YScale','log');
    xlim([0 8]);
    ylim([0 30]);
    legend(hLegend, labelsLegend, ...
           'Location','eastoutside', ...
           'Interpreter','latex', ...
           'FontSize',9, ...
           'Box','on');
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 graphx_cm graphy_cm]);  % 12x12 cm
    filename = sprintf('mu_linear_comparison.png');
    exportgraphics(gcf, filename, 'Resolution', 300, 'BackgroundColor', 'white')
end