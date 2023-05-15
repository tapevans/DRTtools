%% Extract Impedence for DRT from P2D Model Results
clear all; close all; clc;

%% Switch to current directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)
    

%% Parameters
    % FLAG.P2DvsROM = 1; % 1) if from P2D, 2) if from Ho-Kalman
    % filename = 'SS_EIS_SOC50.mat';
    FLAG.P2DvsROM = 0; % 1) if from P2D, 2) if from Ho-Kalman
    filename = 'ROM_SS_SOC50_Ts1.mat';
    foldername = 'Test';
    save_folderpath = [pwd filesep 'myResults' filesep foldername];


%% Load Data
data = load(filename);

%% Convert Results
if FLAG.P2DvsROM
%% Convert Frequencies
    freq           = data.Z_results(:,data.P.SS.omega)/(2*pi);
    Z_prime        = data.Z_results(:,data.P.SS.Z_Re);
    Z_double_prime = data.Z_results(:,data.P.SS.Z_Im);

%% Save Results
    if ~isfolder(save_folderpath)
        mkdir(save_folderpath)
    end
    save([save_folderpath filesep 'INPUTS_' filename], 'freq','Z_prime','Z_double_prime')

else
    %% Define Frequency
        freq_omega = logspace(-1,11,101)/(2*pi);
    for i = 1:length(data.sys_HK)
        for OO = 1:length(data.sys_HK{i})
    
        %% Convert to CT
            sys = d2c(data.sys_HK{i});
    
        %% Get Impedance
            Z = zeros(length(freq_omega),1);
            for k = 1:length(freq_omega)
                s      = freq_omega(k)*(1i);
                Z_temp = sys.C * ((s*eye(size(sys.A)) - sys.A)\sys.B) + sys.D;
                Z(k)   = Z_temp(OO); %%%%%%%% Hard-coded for cell voltage
            end
            A_c = 0.1; %%%% Hardcoded from P2D
            multiple = -A_c^-1;
            Z_new = multiple*Z;
    
    
        %% Get Save Outputs
            freq           = (freq_omega/(2*pi))';
            Z_prime        = real(Z_new);
            Z_double_prime = imag(Z_new);
    
    
        %% Save Results
            if ~isfolder(save_folderpath)
                mkdir(save_folderpath)
            end
            if i < length(data.sys_HK)
                if OO == 1
                    save([save_folderpath filesep 'INPUTS_' data.RESULTS.foldername{i} '_' data.RESULTS.foldername{1} '_' filename], 'freq','Z_prime','Z_double_prime')
                else
                    save([save_folderpath filesep 'INPUTS_' data.RESULTS.foldername{i} '_' data.RESULTS.foldername{i} '_' filename], 'freq','Z_prime','Z_double_prime')
                end
            else
                save([save_folderpath filesep 'INPUTS_AllOutputs' '_' data.RESULTS.foldername{OO} '_' filename], 'freq','Z_prime','Z_double_prime')
            end
        end
    end
end