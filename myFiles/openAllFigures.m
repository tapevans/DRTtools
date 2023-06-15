%% Open All Resulting Figures
clear all; close all; clc;

%% Switch to current directory
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)

%% Parameters
    FLAG.PlotAll = 1;
    FigArrange   = 1;
    FLAG.UsingLaptop = 0; % The issue here was the scale of the monitor was not set to 100%

    foldername = 'Test';
    save_folderpath = [pwd filesep 'myResults' filesep foldername];


%% Plot Pointers
    % Types
    i = 0;
    P.Combined            = i; i = i+1;
    P.Individual_CellVolt = i; i = i+1;
    P.Individual_Desired  = i; i = i+1;

    % Variables
    i = 1;
    P.cell_voltage = i; i = i+1;
    P.delta_phi    = i; i = i+1;
    P.eta          = i; i = i+1;
    P.C_Liion      = i; i = i+1;
    P.C_Li         = i; i = i+1;
    P.delta_C_Li   = i; i = i+1;

    numVar = i - 1;

    DataName{P.cell_voltage} = 'Cell Voltage';
    DataName{P.delta_phi}    = '\Delta\phi';
    DataName{P.eta}          = '\eta';
    DataName{P.C_Liion}      = 'C_{Li^+}';
    DataName{P.C_Li}         = 'C_{Li,surf}';
    DataName{P.delta_C_Li}   = '\Delta C_{Li}';

    ColorVec = [0      0.4470 0.7410
                0.8500 0.3250 0.0980
                0.4940 0.1840 0.5560
                0.4660 0.6740 0.1880
                0.3010 0.7450 0.9330
                0.6350 0.0780 0.1840];


%% Get Filenames of Figures
% Make a list of simulations in the project folder
    oldFolder = cd(save_folderpath);
    list = dir('*.fig*');
    num_files = length(list);
    for j = 1:num_files % Creates a cell array with all simulations' full path name
        if ~exist('sim_filenames')
            sim_filenames{1} = [pwd filesep list(j).name];
        else
            sim_filenames{end+1,1} = [pwd filesep list(j).name];
        end
    end
    %Go back to oldFolder
    cd(oldFolder)


%% Duplicate the Individual Cell Voltage Plot
    for j = 1:num_files
        idx = strfind(sim_filenames{j},'CellVoltage_CellVoltage');
        if ~isempty(idx)
            sim_filenames{end+1,1} = sim_filenames{j};
            list(end+1,1) = list(j,1);
        end
    end
    num_files = num_files+1;


%% Get Titles
% Split Full order model (P2D) and ROM from Ho-Kalman
    FOM_IDX = [];
    ROM_IDX = [];
    for j = 1:num_files
        TitleLabels(j).orig = list(j).name;
        TitleLabels(j).temp = list(j).name(5:end);
        TitleLabels(j).sortIDX = [];
        % Model Type
            if TitleLabels(j).temp(1) == 'S' % SS
                TitleLabels(j).sortIDX = num_files;
                TitleLabels(j).ModelType = 'P2D';
                FOM_IDX = [FOM_IDX,j];
            else
                TitleLabels(j).ModelType = 'ROM: ';
                ROM_IDX = [ROM_IDX,j];
            end
    end

% Split ROM between the combined and individual
    COMB_IDX = [];
    INDV_IDX = [];
    for j = ROM_IDX
        % ROM: Combined or Individual
            if TitleLabels(j).temp(1) == 'A' % AllOutputs
                TitleLabels(j).CorI = 'Combined: ';
                COMB_IDX = [COMB_IDX , j];
            else
                TitleLabels(j).CorI = 'Individual: ';
                INDV_IDX = [INDV_IDX , j];
            end
    end

% For all combined ROM, get title
    for j = COMB_IDX
        idx_first = strfind(TitleLabels(j).temp ,'_');
        idx_ROM   = strfind(TitleLabels(j).temp ,'_ROM');
        varName   = TitleLabels(j).temp(idx_first(1)+1:idx_ROM-1);
        offset    = P.Combined * numVar;
            if strcmp(varName,'C_Li_surf')
                varName = 'C_{Li,surf}';
                TitleLabels(j).sortIDX = offset + P.C_Li;
            elseif strcmp(varName,'C_Liion')
                varName = 'C_{Li^+}';
                TitleLabels(j).sortIDX = offset + P.C_Liion;
            elseif strcmp(varName,'DeltaC_Li')
                varName = '\Delta C_{Li}';
                TitleLabels(j).sortIDX = offset + P.delta_C_Li;
            elseif strcmp(varName,'DeltaPhi')
                varName = '\Delta\phi';
                TitleLabels(j).sortIDX = offset + P.delta_phi;
            elseif strcmp(varName,'eta')
                varName = '\eta';
                TitleLabels(j).sortIDX = offset + P.eta;
            elseif strcmp(varName,'CellVoltage')
                varName = 'Cell Voltage';
                TitleLabels(j).sortIDX = offset + P.cell_voltage;
            else
                disp(['COMB issue with ' num2str(j)])
            end
        % disp(varName)
        TitleLabels(j).VarName = varName;
    end

% For all individual ROM, get title
    for j = INDV_IDX
        idx_ROM      = strfind(TitleLabels(j).temp ,'_ROM');
        temp2        = TitleLabels(j).temp(1:idx_ROM-1);
        idx_cellVolt = strfind(temp2,'CellVoltage');
    % Individual ROM with cell voltage as output
        if ~isempty(idx_cellVolt)
            offset  = numVar * P.Individual_CellVolt;
            varName = TitleLabels(j).temp(1:idx_cellVolt(end)-2);
                if strcmp(varName,'C_Li_surf')
                    varName = 'C_{Li,surf}';
                    TitleLabels(j).sortIDX =  offset + P.C_Li;
                elseif strcmp(varName,'C_Liion')
                    varName = 'C_{Li^+}';
                    TitleLabels(j).sortIDX = offset + P.C_Liion;
                elseif strcmp(varName,'DeltaC_Li')
                    varName = '\Delta C_{Li}';
                    TitleLabels(j).sortIDX = offset + P.delta_C_Li;
                elseif strcmp(varName,'DeltaPhi')
                    varName = '\Delta\phi';
                    TitleLabels(j).sortIDX = offset + P.delta_phi;
                elseif strcmp(varName,'eta')
                    varName = '\eta';
                    TitleLabels(j).sortIDX = offset + P.eta;
                elseif strcmp(varName,'CellVoltage')
                    varName = 'Cell Voltage';
                    TitleLabels(j).sortIDX = offset + P.cell_voltage;
                else
                    disp(['INDV CVolt issue with ' num2str(j)])
                end
            TitleLabels(j).VarName = [varName ': '];
            TitleLabels(j).VarName2 = 'Cell Voltage';
    % Individual ROM desired variable as output
        else
            offset = P.Individual_Desired * numVar;
            varLength = (length(temp2) - 1)/2;
            varName = temp2(1:varLength);
                if strcmp(varName,'C_Li_surf')
                    varName = 'C_{Li,surf}';
                    TitleLabels(j).sortIDX = offset + P.C_Li;
                elseif strcmp(varName,'C_Liion')
                    varName = 'C_{Li^+}';
                    TitleLabels(j).sortIDX = offset + P.C_Liion;
                elseif strcmp(varName,'DeltaC_Li')
                    varName = '\Delta C_{Li}';
                    TitleLabels(j).sortIDX = offset + P.delta_C_Li;
                elseif strcmp(varName,'DeltaPhi')
                    varName = '\Delta\phi';
                    TitleLabels(j).sortIDX = offset + P.delta_phi;
                elseif strcmp(varName,'eta')
                    varName = '\eta';
                    TitleLabels(j).sortIDX = offset + P.eta;
                elseif strcmp(varName,'CellVoltage')
                    varName = 'Cell Voltage';
                    TitleLabels(j).sortIDX = offset + P.cell_voltage;
                else
                    disp(['issue with ' num2str(j)])
                end
            TitleLabels(j).VarName  = [varName ': '];
            TitleLabels(j).VarName2 = varName;
        end
    end

    % Create Plot Titles
        for j = 1:num_files
            TitleLabels(j).Title = [TitleLabels(j).ModelType ...
                                    TitleLabels(j).CorI  ...
                                    TitleLabels(j).VarName  ...
                                    TitleLabels(j).VarName2 ];
        end

    % Fix Duplicate Index
        possible_idx = [P.Individual_Desired * numVar , numVar * P.Individual_CellVolt] + P.cell_voltage;
        if TitleLabels(num_files).sortIDX == possible_idx(1)
            TitleLabels(num_files).sortIDX = possible_idx(2);
        else
            TitleLabels(num_files).sortIDX = possible_idx(1);
        end


%% Sort Filenames
    for j = 1:num_files
        sim_filenames_temp{TitleLabels(j).sortIDX,1} = sim_filenames{j};
    end
    sim_filenames = sim_filenames_temp;

%% Plots
    if FLAG.PlotAll
        %% Open Figures
        for i = 1:num_files
            openfig(sim_filenames{i});
            pause(2.0) % Need this so figure is open before ax is called
            % Save Data
                ax = gca;
                PlotData(i).XData = ax.Children.XData;
                PlotData(i).YData = ax.Children.YData;
                PlotData(i).XLabel = ax.XLabel.String;
                PlotData(i).YLabel = ax.YLabel.String;
            
                for j = 1:num_files
                    if TitleLabels(j).sortIDX == i
                        idx = j;
                    end
                end
                PlotData(i).Title  = TitleLabels(idx).Title;
        end

        %% Modify Figures
        pause(1)
        for i = 1:num_files
            f = figure(i);
            ax = gca;
            % Turn off legend
                ax.Legend.Visible = 'off';
            % Limit x-axis
                if i ~= num_files
                vec = ax.XLim;
                vec(1) = 1e-2;
                ax.XLim = vec;
                end
            % Change Units
                f.Units = 'pixels';
            % Font Size
                ax.FontSize = 12;
            % Add Title
                idx = [];
                for j = 1:num_files
                    if TitleLabels(j).sortIDX == i
                        idx = [idx, j];
                    end
                end
                ax.Title.String = TitleLabels(idx(1)).Title;
            % Change Size
                if ~FigArrange
                    f.Position = [10.6000 363.4000 544 410];
                end

        end

%% Plot Overlapping Figures
% Combined ROM (loglog)
    offset = P.Combined * numVar;
    idx_vec = (1:numVar) + offset;
    figure
    loglog(PlotData(idx_vec(1)).XData , PlotData(idx_vec(1)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(1,:) ,'DisplayName' , DataName{1} )
    hold on 
    for k = 2:numVar
        loglog(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(k,:) ,'DisplayName' , DataName{k} )
    end
    lgn = legend;
    lgn.Location = 'northwest';
    ax = gca;
    ax.XLabel.String = PlotData(i).XLabel;
    ax.YLabel.String = PlotData(i).YLabel;
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.String = 'Combined ROM (loglog)';
    ax.XLim = [1e-2 316];
    % ax.YLim = [1e-20 1];
    
% Combined ROM (semilogx)
    offset = P.Combined * numVar;
    idx_vec = (1:numVar) + offset;
    figure
    semilogx(PlotData(idx_vec(1)).XData , PlotData(idx_vec(1)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(1,:) ,'DisplayName' , DataName{1} )
    hold on 
    for k = 2:numVar
        semilogx(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(k,:) ,'DisplayName' , DataName{k} )
    end
    lgn = legend;
    lgn.Location = 'northwest';
    ax = gca;
    ax.XLabel.String = PlotData(i).XLabel;
    ax.YLabel.String = PlotData(i).YLabel;
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.String = 'Combined ROM (semilogx)';
    ax.XLim = [1e-2 316];
    % ax.YLim = [1e-20 1];

% Combined ROM (semilogx) (Zoomed)
    offset = P.Combined * numVar;
    idx_vec = (1:numVar) + offset;
    figure
    semilogx(PlotData(idx_vec(1)).XData , PlotData(idx_vec(1)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(1,:) ,'DisplayName' , DataName{1} )
    hold on 
    for k = 2:numVar
        semilogx(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(k,:) ,'DisplayName' , DataName{k} )
    end
    lgn = legend;
    lgn.Location = 'northwest';
    ax = gca;
    ax.XLabel.String = PlotData(idx_vec(1)).XLabel;
    ax.YLabel.String = PlotData(idx_vec(1)).YLabel;
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.String = 'Combined ROM (semilogx, zoomed)';
    ax.XLim = [1e-2 316];
    ax.YLim = [0 0.05];

%% Individual Cell Voltages
% Individual ROM Cell Voltage (loglog)
    offset = P.Individual_CellVolt * numVar;
    idx_vec = (1:numVar) + offset;
    figure
    k = 1;
    loglog(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(1,:) ,'DisplayName' , DataName{1} )
    hold on 
    for k = 2:numVar
        loglog(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(k,:) ,'DisplayName' , DataName{k} )
    end
    lgn = legend;
    lgn.Location = 'southeast';
    ax = gca;
    ax.XLabel.String = PlotData(idx_vec(1)).XLabel;
    ax.YLabel.String = PlotData(idx_vec(1)).YLabel;
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.String = 'Individual ROM Cell Voltage';
    ax.XLim = [1e-2 316];
    % ax.YLim = [1e-20 1];

%% Individual Desired Variables
% Individual ROM Cell Voltage (loglog)
    offset = P.Individual_Desired * numVar;
    idx_vec = (1:numVar) + offset;
    figure
    k = 1;
    loglog(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(1,:) ,'DisplayName' , DataName{1} )
    hold on 
    for k = 2:numVar
        loglog(PlotData(idx_vec(k)).XData , PlotData(idx_vec(k)).YData , 'LineWidth' , 2 , 'Color' , ColorVec(k,:) ,'DisplayName' , DataName{k} )
    end
    lgn = legend;
    lgn.Location = 'northwest';
    ax = gca;
    ax.XLabel.String = PlotData(idx_vec(1)).XLabel;
    ax.YLabel.String = PlotData(idx_vec(1)).YLabel;
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.String = 'Individual ROM Desired Variables';
    ax.XLim = [1e-2 316];
    % ax.YLim = [1e-20 1];

%% Arrange Figures
        if FigArrange == 1
            fig = gcf;
            NumFig = fig.Number;
        
            Ncol = 3;
        
            for i = 1:NumFig
                f = figure(i);
                k = mod(i-1,Ncol);
                row = mod(fix((i-1)/Ncol),2);
                if row == 0
                    if FLAG.UsingLaptop
                        r = 930;
                    else
                        r = 575;
                        % r = 540;
                    end
                elseif row == 1
                    if FLAG.UsingLaptop
                        r = 405;
                    else
                        r = 62;
                    end
                end
                if FLAG.UsingLaptop
                    f.Position = [k*575+2000 r 560 420];
                else
                    f.Position = [k*575+15 r 560 420];
                end
            end
        end
    end