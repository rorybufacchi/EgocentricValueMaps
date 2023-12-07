function plotRowAveragesWithSpread_alt(data, spreadWidth, plotTitle, xLabel, yLabel)
    % Function to plot row averages and overlay individual data points
    % with a spread indicating the count of discrete values
    %
    % Inputs:
    %   data       - The input matrix
    %   spreadWidth - Width for spreading points (optional, default = 0.05)
    %   plotTitle  - Title of the plot (optional)
    %   xLabel     - X-axis label (optional)
    %   yLabel     - Y-axis label (optional)

    % Set default values for optional parameters
    if nargin < 2 || isempty(spreadWidth)
% % %         spreadWidth = 0.05;
        spreadWidth = 0.075;
    end
    if nargin < 3
        plotTitle = 'Row Averages and Grouped Data Points';
    end
    if nargin < 4
        xLabel = 'Row Number';
    end
    if nargin < 5
        yLabel = 'Value';
    end

    % Calculate row averages
    rowAverages = mean(data, 2);

    % Create bar graph of row averages
% % %     figure;
% % %     bar(rowAverages);
    hold on;

    % Overlay individual data points with systematic spreading
    [nRows, nCols] = size(data);
    for i = 1:nRows
        % Find unique values and their counts in the row
        [uniqueVals, ~, idx] = unique(data(i, :));
        counts = accumarray(idx(:), 1);

        for j = 1:length(uniqueVals)
            val = uniqueVals(j);
            numPoints = counts(j);
            points = find(data(i, :) == val);

            % Spread out points within the group
            for k = 1:numPoints
                % Systematic offset calculation
                offset = (k - (numPoints + 1) / 2) * spreadWidth;
                x = i + offset;

                % Plot the data point
% % %                 plot(x, data(i, points(k)), 'o');
                scatter(x, data(i, points(k)),200, '.k');
            end
        end
    end

    % Adjust the appearance of the plot
    xlabel(xLabel);
    ylabel(yLabel);
    title(plotTitle);

    % Show the plot
    hold off;
end
