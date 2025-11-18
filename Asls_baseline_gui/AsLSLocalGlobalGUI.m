function AsLSLocalGlobalGUI
% ASLSLocalGlobalGUI
%
% A GUI for Asymmetric Least-Squares (ASLS) baseline correction with 
% both global and local corrections. The baseline is computed by:
%   1) Estimating a global baseline using ASLS.
%   2) Computing local corrections on specified intervals.
%   3) Merging and blending these corrections to form a single smooth baseline.
%
% Steps:
%   1) Load a [nRows x nCols] numeric matrix.
%   2) Specify global ASLS parameters.
%   3) Optionally add local intervals with their own (Lambda, p) parameters.
%   4) Click "Preview Baseline" to view the merged (global+local) baseline 
%      over the mean spectrum.
%   5) (Optional) Enable final smoothing with p = 0.5 and a smoothing Lambda.
%   6) Click "Apply Correction" to baseline-correct all rows in parallel.
%   7) "Export" sends corrected data to the workspace.
%   8) "Export Settings" sends the chosen parameters to the workspace.
%
% Two helper functions below:
%   - Aslsbase_single: improved single-interval ASLS.
%   - GlobalLocalBaseline: merges global and local baselines with smooth transitions.

    bgColor = [1 1 1];
    mainFig = uifigure('Name','ASLS Local + Global Correction',...
        'Color', bgColor, 'Position',[100 100 1000 700]);

    % --- We now use a 2x2 grid instead of 2x3
    mainGrid = uigridlayout(mainFig,[2 2], ...
        'RowHeight',{'1x','fit'}, 'ColumnWidth',{'1.5x','1.5x'}, ...
        'Padding',[10 10 10 10], 'RowSpacing',5, 'ColumnSpacing',5);

    %% (1) AXES: Original Data (Mean)
    previewAxes = uiaxes(mainGrid,'BackgroundColor',bgColor);
    previewAxes.Layout.Row    = 1;
    previewAxes.Layout.Column = 1;
    previewAxes.XLabel.String = 'Channels';
    previewAxes.YLabel.String = 'Intensity';
    previewAxes.Title.String  = 'Original Data (Mean) [size image]';
    grid(previewAxes,'on');

    %% (2) AXES: Corrected Data (Mean)
    correctedAxes = uiaxes(mainGrid,'BackgroundColor',bgColor);
    correctedAxes.Layout.Row    = 1;
    correctedAxes.Layout.Column = 2;
    correctedAxes.XLabel.String = 'Channels';
    correctedAxes.YLabel.String = 'Intensity';
    correctedAxes.Title.String  = 'Corrected Data (Mean)';
    grid(correctedAxes,'on');

    %% Bottom panel for controls
    bottomPanel = uipanel(mainGrid,'BackgroundColor',bgColor);
    bottomPanel.Layout.Row    = 2;
    bottomPanel.Layout.Column = [1 2];

    bottomGrid = uigridlayout(bottomPanel,[1 5]);
    bottomGrid.ColumnWidth   = {'fit','fit','2x','fit','fit'};
    bottomGrid.RowHeight     = {'fit'};
    bottomGrid.Padding       = [5 5 5 5];

    %% (3) Data Panel
    dataPanel = uipanel(bottomGrid,'Title','Data','BackgroundColor',bgColor);
    dataPanel.Layout.Column = 1;
    dataPanel.Layout.Row    = 1;

    dataGrid = uigridlayout(dataPanel,[2 1]);
    dataGrid.RowHeight    = {'fit','fit'};
    dataGrid.ColumnWidth  = {'1x'};

    loadButton = uibutton(dataGrid,'Text','Load Data');
    loadButton.Layout.Row    = 1;
    loadButton.Layout.Column = 1;

    previewButton = uibutton(dataGrid,'Text','Preview Baseline');
    previewButton.Layout.Row    = 2;
    previewButton.Layout.Column = 1;

    %% (4) Global Param Panel
    globalParamPanel = uipanel(bottomGrid,'Title','Global ASLS Params','BackgroundColor',bgColor);
    globalParamPanel.Layout.Column = 2;
    globalParamPanel.Layout.Row    = 1;

    gpGrid = uigridlayout(globalParamPanel,[2 2],...
        'ColumnWidth',{'fit','fit'}, 'RowHeight',{'fit','fit'});

    uilabel(gpGrid,'Text','Lambda','HorizontalAlignment','right');
    globalLambdaField = uieditfield(gpGrid,'numeric','Value',1e6,'Limits',[0 Inf]);

    uilabel(gpGrid,'Text','p','HorizontalAlignment','right');
    globalPField = uieditfield(gpGrid,'numeric','Value',0.001,'Limits',[0 1]);

    %% (5) Interval Table Panel
    intervalPanel = uipanel(bottomGrid,'Title','Local Intervals','BackgroundColor',bgColor);
    intervalPanel.Layout.Column = 3;
    intervalPanel.Layout.Row    = 1;

    intervalGrid = uigridlayout(intervalPanel,[3 2],...
        'RowHeight',{'fit','1x','fit'}, 'ColumnWidth',{'1x','1x'});

    infoLabel = uilabel(intervalGrid,'Text',...
        'Add or remove intervals (Start, End, Lambda, p). Overlaps are averaged.',...
        'WordWrap','on');
    infoLabel.Layout.Row    = 1;
    infoLabel.Layout.Column = 1;

    removeIntervalButton = uibutton(intervalGrid,'Text','Remove Interval');
    removeIntervalButton.Layout.Row    = 1;
    removeIntervalButton.Layout.Column = 2;

    intervalTable = uitable(intervalGrid,'Data',cell(0,4),...
        'ColumnName',{'Start','End','Lambda','p'}, ...
        'ColumnEditable',[true true true true], ...
        'ColumnWidth',{60,60,60,60});
    intervalTable.Layout.Row    = 2;
    intervalTable.Layout.Column = [1 2];
    intervalTable.CellEditCallback = @onIntervalEdit; 

    addIntervalButton = uibutton(intervalGrid,'Text','Add Interval',...
        'ButtonPushedFcn',@onAddInterval);
    addIntervalButton.Layout.Row    = 3;
    addIntervalButton.Layout.Column = [1 2];

    %% (6) Final Smoothing Panel
    smoothPanel = uipanel(bottomGrid,'Title','Final Smoothing','BackgroundColor',bgColor);
    smoothPanel.Layout.Column = 4;
    smoothPanel.Layout.Row    = 1;

    spGrid = uigridlayout(smoothPanel,[2 2],'ColumnWidth',{'fit','fit'},'RowHeight',{'fit','fit'});

    smoothCheck = uicheckbox(spGrid,'Text','Enable Smoothing','Value',false);
    smoothCheck.Layout.Row    = 1;
    smoothCheck.Layout.Column = [1 2];

    smoothLambdaField = uieditfield(spGrid,'numeric','Value',1e5,'Limits',[1 Inf],'Enable','off');
    smoothLambdaField.Layout.Row    = 2;
    smoothLambdaField.Layout.Column = 2;

    labelSmoothLambda = uilabel(spGrid,'Text','Lambda','HorizontalAlignment','right');
    labelSmoothLambda.Layout.Row    = 2;
    labelSmoothLambda.Layout.Column = 1;

    %% (7) Actions Panel
    actionPanel = uipanel(bottomGrid,'Title','Actions','BackgroundColor',bgColor);
    actionPanel.Layout.Column = 5;
    actionPanel.Layout.Row    = 1;

    % Expand to 4 rows so we can fit the new "Export Settings" button
    actGrid = uigridlayout(actionPanel,[4 1],...
        'RowHeight',{'fit','fit','fit','fit'}, 'ColumnWidth',{'1x'});

    applyButton              = uibutton(actGrid,'Text','Apply Correction');
    applyButton.Layout.Row   = 1;
    exportButton             = uibutton(actGrid,'Text','Export to Workspace','Enable','off');
    exportButton.Layout.Row  = 2;

    % NEW: Export Settings Button
    exportParamsButton             = uibutton(actGrid,'Text','Export Settings to Workspace');
    exportParamsButton.Layout.Row  = 3;

    statusLabel             = uilabel(actGrid,'Text','Status: Idle');
    statusLabel.Layout.Row  = 4;

    %% INTERNAL VARIABLES
    dataLoaded    = false;
    originalData  = [];
    correctedData = [];
    waveAxis      = [];
    nRows         = 0;
    nCols         = 0;

    %% CALLBACK ASSIGNMENTS
    loadButton.ButtonPushedFcn           = @onLoad;
    previewButton.ButtonPushedFcn        = @onPreview;
    smoothCheck.ValueChangedFcn          = @onSmoothCheck;
    applyButton.ButtonPushedFcn          = @onApply;
    exportButton.ButtonPushedFcn         = @onExport;
    exportParamsButton.ButtonPushedFcn   = @onExportParams; % New callback
    removeIntervalButton.ButtonPushedFcn = @onRemoveInterval;

    %======================================================================
    %  onLoad
    %======================================================================
    function onLoad(~,~)
        vars = evalin('base','whos');
        names = {vars.name};
        if isempty(names)
            uialert(mainFig,'No variables in base workspace.','Load Error');
            return;
        end
        [sel, ok] = listdlg('ListString',names,'SelectionMode','single',...
            'PromptString','Select a numeric matrix:');
        if ~ok, return; end

        varname = names{sel};
        tmpData = evalin('base',varname);
        if ~isnumeric(tmpData) || ~ismatrix(tmpData)
            uialert(mainFig,'Selected variable is not a numeric 2D matrix.','Invalid Data');
            return;
        end
        originalData = double(tmpData);
        [nRows, nCols] = size(originalData);
        waveAxis = 1:nCols;
        dataLoaded = true;

        % Update the title of the left axes to reflect data size
        previewAxes.Title.String = sprintf('Original Data (Mean) [%dx%d]', nRows, nCols);

        % Clear the other axes
        cla(previewAxes,'reset');
        cla(correctedAxes,'reset');
        correctedAxes.Title.String = 'Corrected Data (Mean)';

        exportButton.Enable = 'off';
        statusLabel.Text    = 'Status: Data loaded';
    end

    %======================================================================
    %  onPreview
    %======================================================================
    function onPreview(~,~)
        if ~dataLoaded
            uialert(mainFig,'Load data first.','No Data');
            return;
        end
        statusLabel.Text = 'Status: Building preview...';
        drawnow;

        yMean = mean(originalData,1)';  

        lamG = globalLambdaField.Value;
        pG   = globalPField.Value;

        tblData    = intervalTable.Data;
        nIntervals = size(tblData,1);
        localIntervals = zeros(nIntervals,2);
        localLambdas   = zeros(nIntervals,1);
        localPs        = zeros(nIntervals,1);
        for i = 1:nIntervals
            localIntervals(i,:) = [round(tblData{i,1}), round(tblData{i,2})];
            localLambdas(i)     = tblData{i,3};
            localPs(i)          = tblData{i,4};
        end

        smoothVal = 10;  % boundary blending
        mergedRow = GlobalLocalBaseline(yMean, lamG, pG, ...
                                        localIntervals, localLambdas, localPs, smoothVal);

        if smoothCheck.Value
            lamS = smoothLambdaField.Value;
            mergedRow = Aslsbase_single(mergedRow, lamS, 0.5);
        end

        cla(previewAxes,'reset'); 
        hold(previewAxes,'on');
        plot(previewAxes, waveAxis, yMean, 'b', 'DisplayName','Mean');
        plot(previewAxes, waveAxis, mergedRow, 'r', 'DisplayName','Merged Baseline');
        legend(previewAxes, 'Location','best');
        hold(previewAxes,'off');

        statusLabel.Text = 'Status: Preview updated';
    end

    %======================================================================
    %  onSmoothCheck
    %======================================================================
    function onSmoothCheck(cb,~)
        if cb.Value
            smoothLambdaField.Enable = 'on';
        else
            smoothLambdaField.Enable = 'off';
        end
    end

    %======================================================================
    %  onApply
    %======================================================================
    function onApply(~,~)
        if ~dataLoaded
            uialert(mainFig,'Load data first.','No Data');
            return;
        end
        lamG = globalLambdaField.Value;
        pG   = globalPField.Value;
        statusLabel.Text = 'Status: Running correction...';

        dlg = uiprogressdlg(mainFig, 'Title', 'Applying Correction',...
            'Message','Please wait, running baseline correction...', ...
            'Indeterminate','on');

        try
            tblData    = intervalTable.Data;
            nIntervals = size(tblData,1);
            localIntervals = zeros(nIntervals,2);
            localLambdas   = zeros(nIntervals,1);
            localPs        = zeros(nIntervals,1);
            for i = 1:nIntervals
                localIntervals(i,:) = [round(tblData{i,1}), round(tblData{i,2})];
                localLambdas(i)     = tblData{i,3};
                localPs(i)          = tblData{i,4};
            end

            smoothVal = 10;

            globalBaseline = zeros(nRows, nCols);
            % Attempt parallel for, otherwise fallback
            try
                parfor r = 1:nRows
                    y = originalData(r,:)';
                    z = GlobalLocalBaseline(y, lamG, pG, ...
                                            localIntervals, localLambdas, localPs, smoothVal);
                    if smoothCheck.Value
                        lamS = smoothLambdaField.Value;
                        z = Aslsbase_single(z, lamS, 0.5);
                    end
                    globalBaseline(r,:) = z';
                end
            catch
                for r = 1:nRows
                    y = originalData(r,:)';
                    z = GlobalLocalBaseline(y, lamG, pG, ...
                                            localIntervals, localLambdas, localPs, smoothVal);
                    if smoothCheck.Value
                        lamS = smoothLambdaField.Value;
                        z = Aslsbase_single(z, lamS, 0.5);
                    end
                    globalBaseline(r,:) = z';
                end
            end

            finalBaseline = globalBaseline;
            correctedData = originalData - finalBaseline;

            cla(correctedAxes, 'reset');
            plot(correctedAxes, waveAxis, mean(correctedData,1), 'r', 'LineWidth',1.2);
            correctedAxes.Title.String = 'Corrected Data (Mean)';

            exportButton.Enable = 'on';
            statusLabel.Text = 'Status: Correction done';
            dlg.Message = 'Baseline correction completed!';
            pause(0.5);
            close(dlg);

        catch ME
            close(dlg);
            uialert(mainFig, ['Error: ', ME.message], 'Error');
            statusLabel.Text = 'Status: Error during correction';
        end
    end

    %======================================================================
    %  onExport
    %======================================================================
    function onExport(~,~)
        if isempty(correctedData)
            uialert(mainFig, 'No corrected data to export.', 'Nothing to Export');
            return;
        end
        prompt = {'Enter variable name to export:'};
        dlgtitle = 'Export Variable';
        dims = [1 50];
        definput = {'correctedData'};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer)
            return; 
        end
        varName = answer{1};
        
        assignin('base', varName, correctedData);
        uialert(mainFig, ['Exported as "', varName, '" to workspace.'], 'Export OK');
    end

    %======================================================================
    %  >>> New Callback: onExportParams <<<
    %======================================================================
    function onExportParams(~, ~)
        if ~dataLoaded
            uialert(mainFig, 'No data loaded yet. Load data before exporting settings.', ...
                    'No Data');
            return;
        end

        % Gather global parameters
        paramStruct.GlobalLambda = globalLambdaField.Value;
        paramStruct.GlobalP      = globalPField.Value;

        % Gather local interval parameters
        tblData = intervalTable.Data;
        if ~isempty(tblData)
            paramStruct.LocalIntervals = cell2mat(tblData(:,1:2));
            paramStruct.LocalLambdas   = cell2mat(tblData(:,3));
            paramStruct.LocalPs        = cell2mat(tblData(:,4));
        else
            paramStruct.LocalIntervals = [];
            paramStruct.LocalLambdas   = [];
            paramStruct.LocalPs        = [];
        end

        % Gather smoothing info
        paramStruct.SmoothingEnabled = smoothCheck.Value;
        if smoothCheck.Value
            paramStruct.SmoothingLambda = smoothLambdaField.Value;
        else
            paramStruct.SmoothingLambda = 0;
        end

        % Prompt user for variable name
        prompt = {'Enter variable name to export the parameters:'};
        dlgtitle = 'Export Parameters';
        dims = [1 50];
        definput = {'baselineCorrectionParams'};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer)
            return; 
        end

        % Export
        varName = answer{1};
        assignin('base', varName, paramStruct);
        uialert(mainFig, ['Exported parameters as "', varName, '" to workspace.'], 'Export OK');
    end

    %======================================================================
    %  onAddInterval
    %======================================================================
    function onAddInterval(~,~)
        if ~dataLoaded
            uialert(mainFig, 'Load data first.', 'No Data');
            return;
        end
        row = {1, min(530, nCols), 1e5, 0.001};
        intervalTable.Data = [intervalTable.Data; row];
        statusLabel.Text = 'Status: New interval row added';
    end

    %======================================================================
    %  onRemoveInterval
    %======================================================================
    function onRemoveInterval(~,~)
        if ~dataLoaded
            uialert(mainFig, 'Load data first.', 'No Data');
            return;
        end
        sel = intervalTable.Selection; 
        if isempty(sel)
            uialert(mainFig,'No row selected to remove.','Remove Interval');
            return;
        end
        rowsToRemove = unique(sel(:,1));
        oldData = intervalTable.Data;
        oldData(rowsToRemove,:) = [];
        intervalTable.Data = oldData;

        statusLabel.Text = 'Status: Interval row(s) removed';
        onPreview();
    end

    %======================================================================
    %  onIntervalEdit
    %======================================================================
    function onIntervalEdit(~, ~)
        if dataLoaded
            onPreview();
        end
    end

end

%% ---------------- SUBFUNCTION: Aslsbase_single ----------------
function [z, WW] = Aslsbase_single(y, lambda, p)
    y = y(:);
    m = length(y);
    if m == 0
        z = [];
        WW = [];
        return;
    end
    D = diff(speye(m), 2);
    w = ones(m, 1);
    z = zeros(m, 1);
    WW = [];
    numIter = 10;
    for iter = 1:numIter
        W = spdiags(w, 0, m, m);
        M = W + lambda * (D' * D) + 1e-12 * speye(m);
        C = chol(M, 'lower');
        temp = C \ (w .* y);
        z = C' \ temp;
        w = p * (y > z) + (1 - p) * (y < z);
        WW = [WW, w]; %#ok<AGROW> 
    end
end

%% ---------------- SUBFUNCTION: GlobalLocalBaseline ----------------
function zMerged = GlobalLocalBaseline(y, globalLambda, globalP, ...
                                       localIntervals, localLambdas, localPs, smoothness)
    n = length(y);
    zGlobal = Aslsbase_single(y, globalLambda, globalP);

    nIntervals = size(localIntervals, 1);
    localMat = nan(nIntervals, n);
    for i = 1:nIntervals
        s = round(localIntervals(i,1));
        e = round(localIntervals(i,2));
        if s > e, [s,e] = deal(e,s); end
        s = max(1, s);
        e = min(n, e);
        if e >= s
            ySub = y(s:e) - zGlobal(s:e);
            if isempty(ySub), continue; end
            zLoc = Aslsbase_single(ySub, localLambdas(i), localPs(i));
            localMat(i, s:e) = zLoc';
        end
    end

    localMean = nanmean(localMat, 1);
    localMean = localMean(:);
    zMerged = zGlobal;
    idx = ~isnan(localMean);
    zMerged(idx) = zGlobal(idx) + localMean(idx);

    if nIntervals > 0 && smoothness > 0
        [~, order] = sort(localIntervals(:,1));
        localIntervals = localIntervals(order, :);
        for i = 1:size(localIntervals,1)-1
            idx1_end   = localIntervals(i,2);
            idx2_start = localIntervals(i+1,1);
            blendLen = min([smoothness, ...
                            idx1_end - localIntervals(i,1) + 1, ...
                            localIntervals(i+1,2) - idx2_start + 1]);
            if blendLen > 0
                idx1 = (idx1_end - blendLen + 1):idx1_end;
                idx2 = idx2_start:(idx2_start + blendLen - 1);
                wBlend = linspace(0, 1, blendLen)';
                zMerged(idx1) = (1 - wBlend) .* zMerged(idx1) + wBlend .* zMerged(idx2);
                zMerged(idx2) = (1 - wBlend) .* zMerged(idx1) + wBlend .* zMerged(idx2);
            end
        end
    end
end

