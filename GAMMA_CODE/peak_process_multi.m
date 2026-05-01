gpufunction [peakResults, B_table] = peak_process_multi(elements, materials, thicknesses)
%% PEAK_PROCESS Process gamma spectroscopy data and calculate B(E,x)
% Ben Taylor-Thornton 2026
% 
%  !!!      This program requires MATLAB Paralell Computing Toolbox     !!!
%  !!! This program requires a NVidia CUDA GPU, VRAM of 4GB or greater  !!!
%
% Inputs:
%   elements    : array of elements [ "Eu-152", "Cs-137" ]
%   materials   : array of arrays of materials for each element eg:
%                   { {"Graded-1", "Graded-2", "Graded-3", "W-PLA High Density"}...
%                     {"Graded-1", "W-PLA High Density"} };
%   thicknesses : array of arrays of arrays of thicknesses of each material 
%                   listed for each element in materials array eg:
%                   { ...
%                   { {"1", "1.5", "3","4.5"}, {"1", "1.5", "2", "3", "4", "7.5"}, {"3"} }, ...
%                   { {"1", "3", "4"} }, ...
%                   };
% Peak Data should be in INTERSPEC CSV FORMAT in folder "peak_data" in the
% same folder as this file, folder and files must be added to MATLAB PATH.
%
% File naming:
% peaks_$Element$ $Thickness$mm $Material$.CSV (case sensitive)
% e.g: peaks_Cs-137 3mm Cu-PLA.CSV
% Stacked Measurements of Same materials should have the TOTAL thickness in
% the filename and have stacking properties defined in stackRules
% containerMap after function definitions.
%
%
% For Graded Shields:
%   Use the gradedLibrary function found after function definitions
%   eg:
%       gradedLibrary("Graded-1") = struct( ...
%       'materials', ["W-PLA High Density", "Cu-PLA"], ...
%       'thicknesses_mm', ["3","1"]);
% File naming:
% peaks_$Element$ n($Thickness$mm $Material$).CSV (case sensitive)
% with n many materials
% e.g: peaks_Eu-152 3mm W-PLA High Density 1mm Cu-PLA.CSV
%
% Puck_Rows
% This program uses "Puck Testing Sheet.csv" for information on shielding
% mechanical properties (Mass, diameter, thickness) and propagates them,
% the format of this file is show below: With D = Diameter, X = Thickness &
% M = Mass
%
% Puck ,Ideal Thickness (cm),Ideal Diameter (cm),R_D,D_1,D_2,D_3,D_mean,R_x,
% x_1,x_2,x_3,x_Mean,R_m,m_1,m_2,m_3,m_Mean
%
% Returns:
%   peakResults : table with processed ln(A/A0), mu, etc.
%   B_table     : table with Monte Carlo propagated B(E,x) <- GPU Required
%
% If no input is given a default is used, defined below.



    %% --- Default Input ---
    if nargin < 1 || isempty(elements)
        elements = ["Eu-152"];
    end
    
    if nargin < 2 || isempty(materials)
        materials = { {"W-PLA High Density"} };
    end
    
    if nargin < 3 || isempty(thicknesses)
        thicknesses = { ...
            { {"1", "3"} }...
        };
    end
    
    %% --- Function Definitions 
    function T = readTwoLineHeaderCSV(filename)
    %% readTwoLineHeaderCSV - Reads a CSV with detected import options
    %
    % Inputs:
    %   filename  - CSV filename to read
    %
    % Output:
    %   T - CSV table data

        % Detect import options automatically
        opts = detectImportOptions(filename);
    
        % Data starts from row 3
        opts.DataLines = [3 Inf];
    
        % Read first two rows as raw cell arrays
        header1 = readcell(filename, 'Range', '1:1');
        header2 = readcell(filename, 'Range', '2:2');
    
        % Ensure both header rows have same width
        nCols = max(numel(header1), numel(header2));
        header1(end+1:nCols) = {''};
        header2(end+1:nCols) = {''};
    
        % Convert to strings
        header1 = string(header1);
        header2 = string(header2);
    
        % Replace missing entries with empty string
        header1(ismissing(header1)) = "";
        header2(ismissing(header2)) = "";
    
        % Combine headers intelligently
        combinedHeaders = strings(1,nCols);
        for i = 1:nCols
            if header1(i) == "" && header2(i) == ""
                combinedHeaders(i) = "Var" + i;
            elseif header1(i) == ""
                combinedHeaders(i) = header2(i);
            elseif header2(i) == ""
                combinedHeaders(i) = header1(i);
            else
                combinedHeaders(i) = header1(i) + "_" + header2(i);
            end
        end
    
        % Make valid + unique MATLAB variable names
        combinedHeaders = matlab.lang.makeValidName(combinedHeaders);
        combinedHeaders = matlab.lang.makeUniqueStrings(combinedHeaders);
    
        % Apply names
        opts.VariableNames = combinedHeaders;
    
        % Read table
        T = readtable(filename, opts);
    
    end
    
    function peak_data = appendPeakDataMulti(peak_data, element, material, thickness, filename)
    %% appendPeakDataMulti - Reads a CSV file with multiple rows and appends them to peak_data
    %
    % Inputs:
    %   peak_data - existing table with 9 columns
    %   element   - string, element name
    %   material  - string, material name (or "BASELINE")
    %   thickness - string, thickness (or "BASELINE")
    %   filename  - CSV filename to read
    %
    % Output:
    %   peak_data - table with new rows appended    
    
        % Read CSV
        data = readTwoLineHeaderCSV(filename);
    
        % Ensure numeric columns are doubles (in case CSV reader returns arrays)
        liveTime    = double(data.LiveTime);
        centroid    = double(data.Centroid_keV);
        photoPeak    = double(data.Photopeak_Energy_keV);
        netArea     = double(data.Net_Area_Counts);
        uncertainty = double(data.Net_Area_Uncertainty);
        fwhm        = double(data.FWHM_keV);
        chi         = double(data.Reduced_Chi_Sqr);
    
        % Number of rows in CSV
        nRows = height(data);
    
        % Create vectors of strings for element/material/thickness
        elements_col   = repmat(string(element), nRows, 1);
        material_col   = repmat(string(material), nRows, 1);
        thickness_col  = repmat(string(thickness), nRows, 1);
        mc_col  = repmat(logical(false), nRows, 1);

    
        % Build table for all rows
        newRows = table( ...
            elements_col, ...
            material_col, ...
            thickness_col, ...
            liveTime, ...
            centroid, ...
            photoPeak, ...
            netArea, ...
            uncertainty, ...
            fwhm, ...
            chi, ...
            mc_col, ...
            'VariableNames', peak_data.Properties.VariableNames ...
        );
    
        % Append to peak_data
        peak_data = [peak_data; newRows];
    end
    
    function puckRows = getPuckRows(puck_data_raw, stackRules, material, thickness_mm)
    %% getPuckRows - Gets relevant puck data for peak_data table
    %
    % Inputs:
    %   puck_data_raw - CSV table data of puck data
    %   stackRules - containerMap, definitions of stacking profiles
    %   material  - string, material name (or "BASELINE")
    %   thickness_mm - string, thickness (or "BASELINE")
    %
    % Output:
    %   puckRows - table with puck data for given material



        puckRows = table();   % Empty return table
        currMat = material;   % Assignment
        
        % Detect if Graded Shield or not
        if startsWith(currMat,"Graded") % Graded Shields
            fprintf("Start Graded Material Import: %s\n", currMat)
            shieldDef = gradedLibrary(currMat);

            % Iterate through Each Graded material from library
            for g = 1:length(shieldDef.materials)
                target_cm = str2double(shieldDef.thicknesses_mm(g)) / 10;
                
                % Try to match to Non-Stacked Graded Materials
                directIdx = strcmp(puck_data_raw.Puck, shieldDef.materials(g)) & ...
                            abs(puck_data_raw.IdealThickness_cm_ - target_cm) < 1e-6;
                
                if any(directIdx)
                    fprintf("Graded, Non-Stacked\n");

                    puckRows = puck_data_raw(directIdx, :);
                    return
                
                
                % Non-Stacked Graded Materials
                else
    
                    key = shieldDef.materials(g) + "|" + string(shieldDef.thicknesses_mm(g));
                    fprintf("Graded, Stack: %s\n", key);
                    if isKey(stackRules, key)
                
                        idealStack = stackRules(key);
                
                        for k = 1:length(idealStack)
                
                            idx = strcmp(puck_data_raw.Puck, material) & ...
                                  abs(puck_data_raw.IdealThickness_cm_ - idealStack(k)) < 1e-6;
                
                            puckRows = [puckRows; puck_data_raw(idx, :)];
                        end
                
                    else
                        error("No puck or stack rule found for %s %s mm",...
                            string(shieldDef.materials(g)), string(shieldDef.thicknesses_mm(g)))
                    end
                end
            end
            

            
        else % Non-Graded Shields
            target_cm = thickness_mm / 10;
            
            % Try direct match
            directIdx = strcmp(puck_data_raw.Puck, material) & ...
                    abs(puck_data_raw.IdealThickness_cm_ - target_cm) < 1e-6;
    
            if any(directIdx)
                puckRows = puck_data_raw(directIdx, :);
                return
            end

            % Try stack rule
            key = material + "|" + string(thickness_mm);
            fprintf(material + "|" + string(thickness_mm))
            fprintf("Non-Graded, Stack: %s\n", key);
            if isKey(stackRules, key)
        
                idealStack = stackRules(key);
        
                for k = 1:length(idealStack)
        
                    idx = strcmp(puck_data_raw.Puck, material) & ...
                          abs(puck_data_raw.IdealThickness_cm_ - idealStack(k)) < 1e-6;
        
                    puckRows = [puckRows; puck_data_raw(idx, :)];
                end
        
            else
                error("No puck or stack rule found for %s %d mm", material, thickness_mm)
            end
        end
    end
    
    function B_table = compute_B_table(peakResults, Nsim)
    %% B_table - Calculates broad-beam factor for a peak data table
    %
    % Inputs:
    %   peakResults - table of peak data
    %   Nsim - Number of Monte-Carlo Simulations to run
    %
    % Output:
    %   B_table - peak data with broad beam factor & uncertainty

        if nargin < 2
            Nsim = 15000; % default MC iterations
        end
    
        % Unique materials and energies
        materials_all = unique(peakResults.Material);
        energies_all  = unique(peakResults.PhotoPeak);
    
        % Initialize B_table with 8 columns
        B_table = table('Size',[0 6], ...
                        'VariableTypes',{'string','double','double',...
                                        'double','double','double'}, ...
                        'VariableNames',{'Material','Thickness','Ideal_Thickness',...
                                        'Energy','B_mean','B_delta'});
    
        % Loop over materials
        for mIdx = 1:length(materials_all)
            matName = materials_all{mIdx};
    
            % Skip BASELINE counts for B(E,x)
            if strcmp(matName, "BASELINE")
                fprintf('Skipping BASELINE material\n');
                continue
            end
    
            % Loop over energies
            for eIdx = 1:length(energies_all)
                E = energies_all(eIdx);
    
                % Select data subset for this material and energy
                idx = peakResults.Material == matName & peakResults.PhotoPeak == E;
                data_subset = peakResults(idx,:);
    
                if isempty(data_subset)
                    continue
                end
    
                % Extract numeric variables
                y  = double(data_subset.ln_transmisson(:)');
                dy = double(data_subset.ln_transmisson_unc(:)');
                x  = double(data_subset.x_total(:)');
                dx = double(data_subset.x_total_unc(:)');
                mu = double(data_subset.mu_linear(:)');
                dmu = double(data_subset.mu_unc(:)');
    
                Npts = length(x);
                if Npts < 4
                    fprintf('%s : %.3f keV | Rejected (N = %d thicknesses < 3)\n', ...
                            matName, E, Npts);
                    continue
                end
    
                % Weighted mu_mean and uncertainty
                mu_mean  = sum(mu ./ dmu.^2) / sum(1 ./ dmu.^2);
                dmu_mean = sqrt(1 / sum(1 ./ dmu.^2));
                fprintf('\n%s : %.3f | Mu_mean = %.2f +- %.5f cm^-1\n', ...
                        matName, E, mu_mean, dmu_mean);
                    
                % Parameters
                chunkSize = 2e7;  % safe for 4 GB VRAM with Npts approx 100-200
                nChunks = ceil(Nsim / chunkSize);
                
                % Initialize accumulators on GPU (single precision)
                B_sum_gpu     = zeros(1, Npts, 'single', 'gpuArray');
                B_sq_sum_gpu  = zeros(1, Npts, 'single', 'gpuArray');
                
                fprintf('Starting fully GPU-resident Monte Carlo simulation...\n');
                
                % Monte Carlo in chunks
                for c = 1:nChunks
                    thisChunk = min(chunkSize, Nsim - (c-1)*chunkSize);
                    
                    tic  % start timing for this chunk
                    
                    % Random numbers directly on GPU (single precision)
                    Y_rand   = (single(y) + single(dy) .* ...
                                randn([thisChunk, Npts], 'single', 'gpuArray'));
                    X_rand   = (single(x) + single(dx) .* ...
                                randn([thisChunk, Npts], 'single', 'gpuArray'));
                    mu_rand  = (single(mu_mean) + single(dmu_mean) .* ...
                                randn([thisChunk,1], 'single', 'gpuArray'));
                    
                    % Monte Carlo propagation
                    Bsim = exp(Y_rand + mu_rand .* X_rand);  % single precision on GPU
                    
                    % Accumulate sums on GPU
                    B_sum_gpu    = B_sum_gpu   + sum(Bsim, 1);
                    B_sq_sum_gpu = B_sq_sum_gpu + sum(Bsim.^2, 1);
                    
                    % Timing
                    wait(gpuDevice);   % <-- THIS IS CRITICAL
                    chunkTime = toc;
                    itersPerSec = thisChunk / chunkTime;
                    
                    fprintf('Chunk %d/%d done | %.2f iterations/sec | chunkTime = %.2f s\n', ...
                        c, nChunks, itersPerSec, chunkTime);
                end
                
                % Bring results back to CPU and compute final statistics
                B_sum_cpu    = gather(B_sum_gpu);
                B_sq_sum_cpu = gather(B_sq_sum_gpu);
                
                B_mean  = double(B_sum_cpu) / Nsim;
                B_delta = sqrt(double(B_sq_sum_cpu) / Nsim - B_mean.^2);
                
                fprintf('\n%s : %.3f | Mu_mean = %.2f +- %.5f cm^-1\n', ...
                    matName, E, mu_mean, dmu_mean);
                    
                for i = 1:Npts
                    thickRaw = data_subset.Thickness(i);  % original entry (string or numeric)
                    if isstring(thickRaw) || ischar(thickRaw)
                        if strcmp(thickRaw,"BASELINE")
                            ideal_thick = NaN;
                        else
                            ideal_thick = str2double(thickRaw);
                        end
                    else
                        ideal_thick = thickRaw;
                    end
    
                    % Append row
                    B_table = [B_table; {matName, x(i), ideal_thick, E, B_mean(i), B_delta(i)}];
                end
            end
        end
    end    


    %% -- DATA ASSIGNMENT
    
    % All data in CM except where stated otherwise!
    
    % Peak_Data Table format, made to match INTERSPEC Peak_Data csv format
    peak_data = table('Size',[0 11], ...
              'VariableTypes',{'string','string','string', 'double', 'double',...
                              'double', 'double', 'double', 'double', 'double' 'logical'}, ...
              'VariableNames',{'Element', ...
                               'Material', ...
                               'Thickness', ...
                               'LiveTime', ...
                               'Centroid', ...
                               'PhotoPeak', ...
                               'Peak_Area', ...
                               'Peak_uncertainty', ...
                               'FWHM', ...
                               'Reduced_chi', ...
                               'Monte_Carlo'});
    gradedLibrary = containers.Map;

    gradedLibrary("Graded-1") = struct( ...
        'materials', ["W-PLA High Density", "Cu-PLA"], ...
        'thicknesses_mm', ["3","1"]); % Materials and thicknesses are ORDERED!
    
    gradedLibrary("Graded-2") = struct( ...
        'materials', ["W-PLA High Density","Cu-PLA"], ...
        'thicknesses_mm', ["3","1.5"]);

    gradedLibrary("Graded-3") = struct( ...
        'materials', ["W-PLA High Density","Cu-PLA"], ...
        'thicknesses_mm', ["3.5","1.5"]);
    
    stackRules = containers.Map;
    % Key format: "Material|Thickness_mm" <-- MM!!
    
    stackRules("W-PLA High Density|4") = [0.3 0.1]; % <-- CM IN ARRAY
    stackRules("W-PLA High Density|3.5") = [0.15 0.2];
    stackRules("W-PLA High Density|7.5") = [0.1 0.15 0.2 0.3];

    stackRules("Cu-PLA|4") = [0.3 0.1];
    stackRules("Pb|2.68") = [0.1321 0.1320 ]; % <-- Workaround
    stackRules("Cu-PLA|4.5") = [0.15 0.3];
    
    %% END OF DATA ASSIGNMENT

    puck_data = readtable("Puck Testing Sheet.csv");  
    puck_data.m_sd = NaN(height(puck_data),1);
    puck_data.m_sigma = NaN(height(puck_data),1);
    puck_data.x_sd = NaN(height(puck_data),1);
    puck_data.x_sigma = NaN(height(puck_data),1);
    puck_data.D_sd = NaN(height(puck_data),1);
    puck_data.D_sigma = NaN(height(puck_data),1);
    
    % Calculate puck data means and SD
    for i = 1:height(puck_data)
        if isnan(puck_data.D_1(i))
            puck_data.D_1(i) = puck_data.D_mean(i);
            puck_data.D_2(i) = puck_data.D_mean(i);
            puck_data.D_3(i) = puck_data.D_mean(i);
            puck_data.x_1(i) = puck_data.x_Mean(i);
            puck_data.x_2(i) = puck_data.x_Mean(i);
            puck_data.x_3(i) = puck_data.x_Mean(i);
            puck_data.m_1(i) = puck_data.m_Mean(i);
            puck_data.m_2(i) = puck_data.m_Mean(i);
            puck_data.m_3(i) = puck_data.m_Mean(i);
        end
        puck_data.D_Mean(i) = mean([puck_data.D_1(i) puck_data.D_2(i) puck_data.D_3(i)]);
        puck_data.x_Mean(i) = mean([puck_data.x_1(i) puck_data.x_2(i) puck_data.x_3(i)]);
        puck_data.m_Mean(i) = mean([puck_data.m_1(i) puck_data.m_2(i) puck_data.m_3(i)]);
        puck_data.D_sd(i) = std([puck_data.D_1(i) puck_data.D_2(i) puck_data.D_3(i)]);
        puck_data.x_sd(i) = std([puck_data.x_1(i) puck_data.x_2(i) puck_data.x_3(i)]);
        puck_data.m_sd(i) = std([puck_data.D_1(i) puck_data.D_2(i) puck_data.D_3(i)]);
        puck_data.D_sigma(i) = sqrt((puck_data.D_sd(i)/sqrt(3))^2 + (puck_data.R_D(i)/6)^2);
        puck_data.x_sigma(i) = sqrt((puck_data.x_sd(i)/sqrt(3))^2 + (puck_data.R_x(i)/6)^2);
        puck_data.m_sigma(i) = sqrt((puck_data.m_sd(i)/sqrt(3))^2 + (puck_data.R_m(i)/6)^2);
    end
    
    %% Extract Peak Data for Each Case
    % Enumerate Elements
    for el = 1:length(elements)
        fprintf('Element: %s\n', elements(el));
        % Append baseline first
        filename = fullfile(pwd, "peak_data", sprintf("peaks_%s Baseline.CSV", elements(el)));
        peak_data = appendPeakDataMulti(peak_data, elements(el), "BASELINE", "BASELINE", filename);
        
        % Get current element's materials and thicknesses
        mats = materials{el};           % cell array of strings
        
        

        
        thicks = thicknesses{el};       % cell array of cell arrays of strings
        
        
        for m = 1:length(mats)
            currMat = mats{m};          % string
            currThicks = thicks{m};     % cell array of strings
            
            if startsWith(currMat,"Graded")
                shieldDef = gradedLibrary(currMat);
                filename = sprintf("peaks_%s", elements(el));
                fprintf('Graded Shield: %s\n# of Materials: %d\nMaterals %s\n', ...
                        currMat ,length(shieldDef.materials), strjoin(shieldDef.materials, ', '));

                for g = 1:length(shieldDef.materials)
            
                    mat_g   = shieldDef.materials(g);
                    thick_g = shieldDef.thicknesses_mm(g);
                    filename = filename + ...
                        sprintf(" %smm %s", thick_g, mat_g);
                end
                
                filename = fullfile(pwd, 'peak_data', sprintf('%s.CSV', filename));
                peak_data = appendPeakDataMulti(peak_data, elements(el), currMat, "", filename);

                
            else 
                for t = 1:length(currThicks)
                    currThick = currThicks{t};   % string
                    
                    fprintf('Material: %s, Thickness: %s\n', currMat, currThick);
                    
                    % Build CSV filename
                    filename = fullfile(pwd, 'peak_data', sprintf("peaks_%s %smm %s.CSV", ...
                                                                   elements(el), ...
                                                                   currThick, ...
                                                                   currMat));
                    
                    % Append row(s) to peak_data
                    peak_data = appendPeakDataMulti(peak_data, elements(el), currMat, ...
                                                    currThick, filename);
                end
            end
        end
    end
    
    %% Calculation of mu, mu_mass and more
    peak_data.BaselineArea        = NaN(height(peak_data),1);
    peak_data.BaselineUncertainty = NaN(height(peak_data),1);
    
    for i = 1:height(peak_data)
    
        if peak_data.Material(i) ~= "BASELINE"
    
            idx = peak_data.Element == peak_data.Element(i) & ...
                  peak_data.Material == "BASELINE" & ...
                  peak_data.Thickness == "BASELINE" & ...
                  peak_data.PhotoPeak == peak_data.PhotoPeak(i);
    
            row = peak_data(idx,:);
    
            if ~isempty(row)
                peak_data.BaselineArea(i)        = row.Peak_Area;
                peak_data.BaselineUncertainty(i) = row.Peak_uncertainty;
            end
        end
    end
    
    % Pre-allocate columns
    peak_data.mu_linear = NaN(height(peak_data),1);
    peak_data.mu_unc = NaN(height(peak_data),1);
    peak_data.mu_int = NaN(height(peak_data),1);
    peak_data.mu_mass_unc = NaN(height(peak_data),1);
    peak_data.x_total = NaN(height(peak_data),1);
    peak_data.x_total_unc = NaN(height(peak_data),1);
    peak_data.ln_transmisson = NaN(height(peak_data),1);
    peak_data.ln_transmisson_unc = NaN(height(peak_data),1);
    
    for i = 1:height(peak_data)
    
        % Skip baseline
        if peak_data.Material(i) == "BASELINE"
            continue
        end
    
        % --- Get baseline peak ---
        idxBase = peak_data.Element == peak_data.Element(i) & ...
                  peak_data.Material == "BASELINE" & ...
                  peak_data.Thickness == "BASELINE" & ...
                  peak_data.PhotoPeak == peak_data.PhotoPeak(i);
    
        if ~any(idxBase)
            warning("No baseline found for row %d", i);
            continue
        end
    
        A_baseline = peak_data.Peak_Area(idxBase);
        A_measured = peak_data.Peak_Area(i);
    
        deltaA = peak_data.Peak_uncertainty(i);
        deltaA0 = peak_data.Peak_uncertainty(idxBase);
    
        % --- Get puck rows ---
        thickness_mm = str2double(peak_data.Thickness(i));
        puckRows = getPuckRows(puck_data, stackRules, ...
                                peak_data.Material(i), thickness_mm);
        
        % --- Adds total thickness for Graded Shields
        val = str2double(peak_data.Thickness(i));
        
        if isnan(val)
            matKey = peak_data.Material(i);  % e.g. "Graded-3"
            thicknesses = gradedLibrary(matKey).thicknesses_mm;
            peak_data.Thickness(i) = sum(str2double(thicknesses));
        else
            peak_data.Thickness(i) = val;
        end
        
        % Actual Calculations
        N = height(puckRows); % Number of pucks
        peak_data.ln_transmisson(i) = log(A_measured/A_baseline);
        peak_data.ln_transmisson_unc(i) = sqrt( (peak_data.Peak_uncertainty(i)/A_measured )^2 + ...
                                                (deltaA0/A_baseline)^2 );
        x_total = sum(puckRows.x_Mean);
        x_total_delta = 1/height(puckRows) * sum(puckRows.x_Mean);
        mu_linear = -1/x_total * log(A_measured / A_baseline);
        rho_stack = ( ( 4*sum(puckRows.m_Mean ) / ...
                    (pi * ( (1/N) * sum(puckRows.D_Mean) )^2 * sum(puckRows.x_Mean)) ) );
        mu_mass = mu_linear / rho_stack;
        peak_data.mu_linear(i) = mu_linear;
        peak_data.mu_mass(i) = mu_mass;
        peak_data.x_total(i) = x_total;
        peak_data.x_total_unc(i) = x_total_delta;
        % Propagate uncertainty
        mu_linear_delta = sqrt( 1/x_total^4 * ...
                                ( (sum(puckRows.x_sigma.^2)) * ...
                                  (log(A_baseline/A_measured)^2 + ...
                                  (sum(puckRows.x_Mean)^2 * ...
                                  ( (deltaA/A_measured)^2 + ...
                                  (deltaA0/A_baseline)^2 ) ) ...
                                  )));

        rho_stack_delta = rho_stack * sqrt ( ( sqrt(sum(puckRows.m_sigma.^2) / ...
                                                    sum(puckRows.x_Mean) ) )^2 + ...
                                             ( sqrt(sum(puckRows.x_sigma.^2) / ...
                                                    sum(puckRows.x_Mean) ) )^2 + ...
                                             ( 2*sqrt(sum(puckRows.D_sigma.^2) / ...
                                                    sum(puckRows.D_Mean) ) )^2);
        
        mu_mass_delta = (1/rho_stack) * sqrt( mu_linear_delta^2 + ...
                        ((mu_linear*rho_stack_delta) / rho_stack)^2 );
        
        peak_data.mu_mass_unc(i) = mu_mass_delta;
        peak_data.mu_unc(i) = mu_linear_delta;
    
        fprintf(['\nProcessed Peak:\n' ...
         '%s - %s (%s mm) N=%d\n' ...        % Element - Material (thickness)
         '%0.2f keV\n' ...              % Centroid energy
         'Mu: %0.4f+/-%0.4fcm^-1, Delta_Mu: %0.4f\n' ...   % Linear attenuation
         'Mu_m: %0.4f+/-%0.4fcm^2/g, Delta_Mu_m: %0.4f\n' ... % Mass attenuation
         'Rho_stack: %0.4f+/-%0.4fgcm^-3, Delta_rho: %0.4f \n'], ... % density of puck(s)
         peak_data.Element(i), ...
         peak_data.Material(i), ...
         peak_data.Thickness(i), ...
         N, ...
         peak_data.Centroid(i), ...
         mu_linear, ...
         peak_data.mu_unc(i), ...
         ((mu_linear_delta/mu_linear) * 100), ... 
         mu_mass, ...
         mu_mass_delta, ...
         ((mu_mass_delta/mu_mass) * 100), ...
         rho_stack, ...
         (rho_stack_delta), ...
         ((rho_stack_delta/rho_stack)*100));  
    
    end
    peakResults = peak_data;

    %% Calculation of Broad Beam Factor
    %B_table = compute_B_table(peakResults, 1.5e8);
    B_table = {};
end