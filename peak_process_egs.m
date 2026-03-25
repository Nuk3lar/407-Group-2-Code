function [egs_peakResults] = peak_process_egs()
%% PEAK_PROCESS Process gamma spectroscopy
% Ben Taylor-Thornton 2026
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
%
% Returns:
%   egs_peakResults : table with processed ln(A/A0), mu, etc.
%   
    elements = ["Cs-137", "Eu-152", "Am-241"];
    nsim = [10, 25, 25];

    materials = { {"W"}, {"W"}, {"Cu"} };
    % IN MM
    thicknesses = { ...
        { { "1", "2", "3", "4", "5", "7.5" } },...
        { { "1", "2", "3", "3.5", "4", "5", "7.5" } }, ...
        { { "1", "1.5", "2", "3", "4", "5", "7.5" } }...
    };

    function peak_data = appendPeakDataMulti(peak_data, element, material, thickness, filename)
        % --- User parameters ---
        tolerance = 0.005; % MeV, absolute tolerance for peak matching
    
        % --- Read raw particle data ---
        data_raw = readtable(filename);
        E = data_raw.E_MeV;
    
        % --- Histogram edges (0-1 MeV, 2000 bins) ---
        edges = linspace(0,1,2001);
        bin_centers = (edges(1:end-1) + edges(2:end))/2;
        bin_centers = bin_centers(:);
    
        % --- Counts ---
        counts = histcounts(E, edges);
        counts = counts(:);
    
        % --- Smooth --- 
        counts_smooth = movmean(counts, 5); 
        % --- Peak finding --- 
        [pks, locs_idx] = findpeaks(counts_smooth, ...
            'MinPeakProminence', max(counts)*0.01);
    
        nPeaks = length(locs_idx);
    
        fprintf('Detected %d peaks for %s, %s\n', nPeaks, element, material);
    
        % --- Preallocate ---
        centroid        = zeros(nPeaks,1);
        netArea         = zeros(nPeaks,1);
        uncertainty     = zeros(nPeaks,1);
        matchedElement  = strings(nPeaks,1);
        matchedEnergy   = NaN(nPeaks,1);
    
        window = 5; % bins around peak for centroid/area
    
        for i = 1:nPeaks
            idx0 = locs_idx(i);
            idx = (idx0-window):(idx0+window);
            idx = idx(idx>0 & idx<=length(counts));
    
            if sum(counts(idx)) > 0
                % Weighted centroid (MeV)
                centroid(i) = sum(bin_centers(idx).*counts(idx))/sum(counts(idx))*1000; % keV
                % Peak area
                netArea(i) = sum(counts(idx));
                % Uncertainty (Poisson)
                uncertainty(i) = sqrt(sum(counts(idx)));
            else
                centroid(i) = NaN;
                netArea(i) = 0;
                uncertainty(i) = NaN;
            end
    
            % --- Peak matching ---
            bestMatch = "";
            bestEnergy = NaN;
            bestDiff = inf;
            for j = 1:length(peakLibrary)
                E_lib = peakLibrary(j).Energy; % MeV
                diff = abs(centroid(i)/1000 - E_lib); % convert keV→MeV
                if diff < tolerance && diff < bestDiff
                    bestDiff = diff;
                    bestMatch = peakLibrary(j).Element;
                    bestEnergy = E_lib;
                end
            end
    
            if bestMatch == ""
                matchedElement(i) = "Unknown";
            else
                matchedElement(i) = bestMatch;
                matchedEnergy(i) = bestEnergy*1000; % store in keV
            end
    
            fprintf('Peak %d: Centroid = %.3f keV, Area = %.1f, Match = %s (%.1f keV)\n', ...
                    i, centroid(i), netArea(i), matchedElement(i), matchedEnergy(i));
        end
    
        % --- Fill remaining table columns ---
        liveTime  = NaN(nPeaks,1);
        photoPeak = matchedEnergy;
        fwhm      = NaN(nPeaks,1);
        chi       = NaN(nPeaks,1);
    
        if string(material) == "BASELINE"
            material_str = "BASELINE";
        else
            material_str = "(EGS) " + string(material) + "-PLA";
        end
        elements_col  = matchedElement;
        material_col  = repmat(material_str, nPeaks,1);
        thickness_col = repmat(string(thickness), nPeaks,1);
        mc_col        = repmat(logical(true), nPeaks,1);
    
        % --- Build table ---
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
    
        % --- Append ---
        peak_data = [peak_data; newRows];
    end
    
    

    %% -- DATA ASSIGNMENT
    
    peakLibrary = struct();

    peakLibrary(1).Element = "Cs-137";
    peakLibrary(1).Energy  = 0.66166;   % MeV
    peakLibrary(2).Element = "Eu-152";
    peakLibrary(2).Energy  = 0.1217800;   % MeV

    peakLibrary(3).Element = "Eu-152";
    peakLibrary(3).Energy  = 0.244700;   % MeV

    peakLibrary(4).Element = "Eu-152";
    peakLibrary(4).Energy  = 0.344280;   % MeV

    peakLibrary(5).Element = "Eu-152";
    peakLibrary(5).Energy  = 0.77890;   % MeV

    peakLibrary(6).Element = "Eu-152";
    peakLibrary(6).Energy  = 0.964080;   % MeV

    peakLibrary(7).Element = "Am-241";
    peakLibrary(7).Energy  = 0.05954;
    
    peakLibrary(8).Element = "Am-241";
    peakLibrary(8).Energy  = 0.02635;% MeV
    % All data in CM except where stated otherwise!
    
    % Peak_Data Table format, made to match INTERSPEC Peak_Data csv format
    peak_data = table('Size',[0 11], ...
              'VariableTypes',{'string','string','string', 'double', 'double',...
                              'double', 'double', 'double', 'double', 'double', 'logical'}, ...
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


    
    %% Extract Peak Data for Each Case
    % Enumerate Elements
    for el = 1:length(elements)
        fprintf('Element: %s\n', elements(el));
        % Append baseline first
        filename = fullfile(pwd, 'EGS_C27_SIM', elements(el), 'Baseline_NS', sprintf('%s_NS_%d_Mil.csv', elements(el), nsim(el)));
        peak_data = appendPeakDataMulti(peak_data, elements(el), "BASELINE", "BASELINE", filename);
        
        % Get current element's materials and thicknesses
        mats = materials{el};           % cell array of strings
        
        

        
        thicks = thicknesses{el};       % cell array of cell arrays of strings
        
        
        for m = 1:length(mats)
            currMat = mats{m};          % string
            currThicks = thicks{m};     % cell array of strings
            
            
            for t = 1:length(currThicks)
                currThick = currThicks{t};   % string
                
                fprintf('Material: %s, Thickness: %s\n', currMat, currThick);
                
                % Build CSV filename

                filename = fullfile(pwd, 'EGS_C27_SIM', ...
                                    elements(el), sprintf('%smm_%s', ...
                                    currThick, currMat), ...
                                    sprintf('%d_Mil_Detector_Log_%s_%s.csv', nsim(el), currThick, currMat));

                
                % Append row(s) to peak_data
                peak_data = appendPeakDataMulti(peak_data, elements(el), currMat, ...
                                                currThick, filename);
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
    
            if height(row) > 1
                warning('Multiple baseline rows found for Element=%s, Peak=%.2f', ...
                        peak_data.Element(i), peak_data.PhotoPeak(i));
            end
            
            if ~isempty(row)
                peak_data.BaselineArea(i)        = row.Peak_Area(1);
                peak_data.BaselineUncertainty(i) = row.Peak_uncertainty(1);
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
        thickness_cm = 0.1 * thickness_mm
        
        % Actual Calculations
        peak_data.ln_transmisson(i) = log(A_measured/A_baseline);
        peak_data.ln_transmisson_unc(i) = sqrt( (peak_data.Peak_uncertainty(i)/A_measured )^2 + ...
                                                (deltaA0/A_baseline)^2 );
        
        mu_linear = -1/thickness_cm * log(A_measured / A_baseline);
        rho_stack = 8.6947; % Hardcoded for W-PLA
        mu_mass = mu_linear / rho_stack;
        peak_data.mu_linear(i) = mu_linear;
        peak_data.mu_mass(i) = mu_mass;
        peak_data.x_total(i) = thickness_cm;
        peak_data.x_total_unc(i) = NaN;
        % Propagate uncertainty
        mu_linear_delta = sqrt( (-(deltaA)/(A_measured*thickness_cm))^2 + ...
                                (-(deltaA0)/(A_baseline*thickness_cm))^2);

        rho_stack_delta = 0;
        
        mu_mass_delta = (1/rho_stack) * sqrt( mu_linear_delta^2 + ...
                        ((mu_linear*rho_stack_delta) / rho_stack)^2 );
        
        peak_data.mu_mass_unc(i) = mu_mass_delta;
        peak_data.mu_unc(i) = mu_linear_delta;
    

    end
    peak_data(strcmp(peak_data.Element, "Unknown"), :) = [];
    egs_peakResults = peak_data;
end