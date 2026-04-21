%% Parameters
energy = 120;       % kV for the CSV file
filter = 4;      % number of filters
bins = 2000;
tube_current = 2; % mA
diameter = 5; %cm4
shielding_1 = 0.0;
shielding_2 = 0.6;
filename = sprintf('XStrahl Data//%dkvP_Filter_%d_3.csv', energy, filter);

% Read CSV
T = readtable(filename);
T.Var4 = NaN(height(T),1);

bin_width = energy / bins;
nelec = (tube_current / 1000) / 1.602e-19;
area = (diameter^2 * pi) / 4; 

% Read XCOM data
filename_xcom = "XCOM Data//XCOM_W-PLA.csv";
opts = detectImportOptions(filename_xcom,'NumHeaderLines',6);
xcomTable = readtable(filename_xcom, opts);

% XCOM energies and mass attenuation coefficients
xcom_energy = xcomTable{:,1};   % MeV
mu_xcom     = xcomTable{:,3};   % cm^2/g

density = 8.8;                  % g/cm^3

% Linear attenuation coefficient
mu_linear = mu_xcom .* density; % 1/cm

% --- Sort data ---
[xcom_energy, sortIdx] = sort(xcom_energy);
mu_linear = mu_linear(sortIdx);

% --- Find duplicate energies (absorption edges) ---
dup_idx = find(diff(xcom_energy) == 0);

% --- Define segment boundaries ---
edges = [1; dup_idx+1; length(xcom_energy)+1];

% --- Allocate output ---
T.Var3 = NaN(size(T.Var1));

% --- Piecewise log-log interpolation ---
for i = 1:length(edges)-1
    idx_range = edges(i):(edges(i+1)-1);
    
    x_seg = xcom_energy(idx_range);
    mu_seg = mu_linear(idx_range);
    
    % Ensure uniqueness within segment
    [x_seg, ia] = unique(x_seg, 'stable');
    mu_seg = mu_seg(ia);
    
    % Points in this segment
    mask = T.Var1 >= min(x_seg) & T.Var1 <= max(x_seg);
    
    if any(mask)
        T.Var3(mask) = exp(interp1( ...
            log(x_seg), ...
            log(mu_seg), ...
            log(T.Var1(mask)), ...
            'linear'));
    end
end

% --- Optional: handle extrapolation ---
below = T.Var1 < xcom_energy(1);
above = T.Var1 > xcom_energy(end);

T.Var3(below) = T.Var3(find(~isnan(T.Var3),1,'first'));
T.Var3(above) = T.Var3(find(~isnan(T.Var3),1,'last'));
% Interpolate mu_linear onto energies in T.Var1


% Read XCOM data
filename_xcom = "XCOM Data//XCOM_Cu-PLA.csv";
opts = detectImportOptions(filename_xcom,'NumHeaderLines',6);
xcomTable = readtable(filename_xcom, opts);

% XCOM energies and mass attenuation coefficients
xcom_energy = xcomTable{:,1};   % MeV
mu_xcom     = xcomTable{:,3};   % cm^2/g

density = 4.77;                  % g/cm^3

% Linear attenuation coefficient
mu_linear = mu_xcom .* density; % 1/cm

% --- Sort data ---
[xcom_energy, sortIdx] = sort(xcom_energy);
mu_linear = mu_linear(sortIdx);

% --- Find duplicate energies (absorption edges) ---
dup_idx = find(diff(xcom_energy) == 0);

% --- Define segment boundaries ---
edges = [1; dup_idx+1; length(xcom_energy)+1];

% --- Allocate output ---
T.Var4 = NaN(size(T.Var1));

% --- Piecewise log-log interpolation ---
for i = 1:length(edges)-1
    idx_range = edges(i):(edges(i+1)-1);
    
    x_seg = xcom_energy(idx_range);
    mu_seg = mu_linear(idx_range);
    
    % Ensure uniqueness within segment
    [x_seg, ia] = unique(x_seg, 'stable');
    mu_seg = mu_seg(ia);
    
    % Points in this segment
    mask = T.Var1 >= min(x_seg) & T.Var1 <= max(x_seg);
    
    if any(mask)
        T.Var4(mask) = exp(interp1( ...
            log(x_seg), ...
            log(mu_seg), ...
            log(T.Var1(mask)), ...
            'linear'));
    end
end

% --- Optional: handle extrapolation ---
below = T.Var1 < xcom_energy(1);
above = T.Var1 > xcom_energy(end);

T.Var4(below) = T.Var4(find(~isnan(T.Var3),1,'first'));
T.Var4(above) = T.Var4(find(~isnan(T.Var3),1,'last'));


countspersecond = 0;
countspersecondshielded = 0;
for i = 1:(length(T.Var2)/3)
    y_value = T.Var2(i*3);
    cps_per_bin = y_value * nelec * bin_width * area;
    countspersecond = countspersecond + cps_per_bin;
    shielded_cps_1 = cps_per_bin * exp(-T.Var4(i*3)*shielding_1);
    shielded_cps_2 = shielded_cps_1 * exp(-T.Var3(i*3)*shielding_2);
    countspersecondshielded = countspersecondshielded + shielded_cps_2;
    fprintf("CPS @ %dkeV: %d bq, %d bq (shielded)\n", T.Var1(i*3)*1000, cps_per_bin, shielded_cps_2)
    
end
format short
fprintf("Bin width: %d ev\n", bin_width)

fprintf("Final CPS: %d bq\n", countspersecond)
fprintf("Shielded CPS: %d bq\n", countspersecondshielded)