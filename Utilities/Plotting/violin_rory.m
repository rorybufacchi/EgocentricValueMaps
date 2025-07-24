%__________________________________________________________________________
% violin.m - Simple violin plot using matlab default kernel density estimation
% Last update: 10/2015
%__________________________________________________________________________
% This function creates violin plots based on kernel density estimation
% using ksdensity with default settings. Please be careful when comparing pdfs
% estimated with different bandwidth!
%
% Differently to other boxplot functions, you may specify the x-position.
% This is usefule when overlaying with other data / plots.
%__________________________________________________________________________
%
% Please cite this function as:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
%
%__________________________________________________________________________
%
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
% bw:        Kernel bandwidth. (default []); prescribe if wanted as follows:
%            a) if bw is a single number, bw will be applied to all
%            columns or cells
%            b) if bw is an array of 1xm or mx1, bw(i) will be applied to cell or column (i).
%            c) if bw is empty (default []), the optimal bandwidth for
%            gaussian kernel is used (see Matlab documentation for
%            ksdensity()
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% bw:    bandwidth of kernel
%__________________________________________________________________________
%{
% Example1 (default):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
[h,L,MX,MED]=violin(Y);
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example2 (specify facecolor, edgecolor, xlabel):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
violin(Y,'xlabel',{'a','b','c','d'},'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','r--')
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example3 (specify x axis location):

disp('this example uses the statistical toolbox')
Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
violin(Y,'x',[-1 .7 3.4 8.8],'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none',...
'bw',0.3,'mc','k','medc','r-.')
axis([-2 10 -0.5 20])
ylabel('\Delta [yesno^{-2}]','FontSize',14)

%Example4 (Give data as cells with different n):

disp('this example uses the statistical toolbox')

Y{:,1}=rand(10,1);
Y{:,2}=rand(1000,1);
violin(Y,'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none','bw',0.1,'mc','k','medc','r-.')
ylabel('\Delta [yesno^{-2}]','FontSize',14)
%}
%%
function[h,L,MX,MED,bw]=violin_rory(Y,varargin) % Renamed function

%defaults:
%_____________________
xL=[];
fc=[1 0.5 0];
lc='k';
alp=0.5;
mc='k';
medc='r';
b=[]; %bandwidth
plotlegend=1;
plotmean=1;
plotmedian=1;
x = [];
%_____________________

%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end

%get additional input parameters (varargin)
p = inputParser;
p.KeepUnmatched = true; % Allow other parameters not specified here
addParameter(p, 'xlabel', []);
addParameter(p, 'facecolor', [1 0.5 0]);
addParameter(p, 'edgecolor', 'k');
addParameter(p, 'facealpha', 0.5);
addParameter(p, 'mc', 'k', @(c)ischar(c) || isempty(c) || isnumeric(c)); % Allow [], colorspec, 'none'
addParameter(p, 'medc', 'r', @(c)ischar(c) || isempty(c) || isnumeric(c)); % Allow [], colorspec, 'none'
addParameter(p, 'bw', [], @isnumeric);
addParameter(p, 'plotlegend', true, @islogical);
addParameter(p, 'x', [], @isnumeric);
addParameter(p, 'support', [-Inf, Inf], @(s)isnumeric(s) && numel(s)==2 && s(1)<s(2)); % Added support option
addParameter(p, 'PlotDataPoints', false, @islogical);
addParameter(p, 'DataPointSize', 8, @isscalar);
addParameter(p, 'DataPointColor', [0.3 0.3 0.3], @(c)ischar(c) || (isnumeric(c) && (numel(c)==3 || size(c,2)==3))); % Allow color name or RGB
addParameter(p, 'DataPointMarker', '.', @ischar);
addParameter(p, 'JitterWidth', 0.8, @(x)isscalar(x) && x>=0 && x<=1); % Max jitter relative to violin half-width

parse(p, varargin{:});

plotDataPoints = p.Results.PlotDataPoints;
dataPointSize = p.Results.DataPointSize;
dataPointColor = p.Results.DataPointColor;
dataPointMarker = p.Results.DataPointMarker;
jitterWidth = p.Results.JitterWidth; % Controls how much horizontal space jitter uses
xL = p.Results.xlabel;
fc = p.Results.facecolor;
lc = p.Results.edgecolor;
alp = p.Results.facealpha;
mc = p.Results.mc;
medc = p.Results.medc;
b = p.Results.bw;
plotlegend = p.Results.plotlegend;
x = p.Results.x;
data_support = p.Results.support; % Get the support range

% Determine if mean/median should be plotted
plotmean = ~isempty(mc) && ~isequal(lower(mc), 'none');
plotmedian = ~isempty(medc) && ~isequal(lower(medc), 'none');


if size(fc,1)==1 && size(Y,2) > 1 % Check if Y has multiple columns before repmat
    fc=repmat(fc,size(Y,2),1);
end

if ~isempty(b)
    if isscalar(b)
        disp(['Using same bandwidth bw = ',num2str(b),' for all groups'])
        b=repmat(b,1, size(Y,2)); % Ensure row vector if needed
    elseif numel(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('Please provide only one bandwidth or an array of b with same length as the number of groups in the data set')
    end
end

% Preallocate outputs
nGroups = size(Y,2);
MX = NaN(1, nGroups);
MED = NaN(1, nGroups);
bw = NaN(1, nGroups);
F = []; % Will be sized later
U = []; % Will be sized later
h = gobjects(1, nGroups); % Use gobjects for graphics handles
p = gobjects(1, 2); % Handles for mean/median lines for legend
FiniteData = cell(1, nGroups); % To store data for plotting points

%% Calculate the kernel density for each group
ksdensity_options = {'support', data_support}; % Use specified support
U_rows = 100; % Default/initial guess for number of points

for i=1:nGroups

    %fprintf('Processing Violin %d...\n', i); % Debug trace
    current_data = Y{i};

    % Ensure data is numeric column vector
    if ~isnumeric(current_data) || ~isvector(current_data)
       warning('Violin %d: Data is not a numeric vector. Skipping.', i);
       if i==1, [F, U] = deal(NaN(U_rows, nGroups)); end % Allocate if first
       continue;
    end
    current_data = current_data(:); % Ensure column vector

    % 1. Filter out non-finite values (NaN, Inf)
    finite_idx = isfinite(current_data);
    if ~all(finite_idx)
        num_removed = sum(~finite_idx);
        % fprintf('Violin %d: Found and removed %d non-finite values.\n', i, num_removed);
        current_data = current_data(finite_idx);
    end

    % 2. Check if any data remains after filtering
    if isempty(current_data)
         fprintf('Violin %d: No finite data points remaining after filtering. Skipping ksdensity.\n', i);
         if i==1 && isempty(F), [F, U] = deal(NaN(U_rows, nGroups)); end % Allocate if first
         if ~isempty(F) % If F/U exist, ensure columns match size
             F(:,i) = NaN(size(F,1), 1);
             U(:,i) = NaN(size(U,1), 1);
         end
         MED(i) = nanmedian(Y{i}(:)); % Use original data for stats
         MX(i) = nanmean(Y{i}(:));
         bw(i) = NaN;
         continue; % Skip to the next iteration
    end

    FiniteData{i} = current_data; % Store the finite data for later plotting

    % 3. Clamp data strictly to the specified support range
    lower_bound = data_support(1);
    upper_bound = data_support(2);
    % Initial clamp:
    current_data_clamped = max(lower_bound, min(upper_bound, current_data));

    % --- !!! NEW: Apply tiny epsilon nudge away from exact boundaries !!! ---
    nudge_eps = 1e-10; % A very small number
    current_data_clamped(current_data_clamped <= lower_bound) = lower_bound + nudge_eps;
    current_data_clamped(current_data_clamped >= upper_bound) = upper_bound - nudge_eps;
    % --- End Epsilon Nudge ---

    % Calculate stats on original finite data (before clamping/nudge)
    MED(i)=median(current_data);
    MX(i)=mean(current_data);

    % --- !!! CRITICAL DEBUG CHECK (uses original bounds for check) !!! ---
    min_val = min(current_data_clamped);
    max_val = max(current_data_clamped);
    tolerance = 1e-9; % Tolerance for floating point comparison
    fail_check = false;
    % Check against the ORIGINAL support bounds, but expect nudge to keep it slightly inside
    if min_val < lower_bound - tolerance % Should definitely not happen now
        fprintf('DEBUG Violin %d: Nudged min value is STILL < lower bound (%f): %f\n', i, lower_bound, min_val);
        fail_check = true;
    elseif min_val <= lower_bound % Check if it's still exactly on the boundary after nudge
         fprintf('DEBUG Violin %d: Nudged min value is STILL <= lower bound (%f): %f\n', i, lower_bound, min_val);
         % This might be okay, but unexpected after nudge
    end
     if max_val > upper_bound + tolerance % Should definitely not happen now
        fprintf('DEBUG Violin %d: Nudged max value is STILL > upper bound (%f): %f\n', i, upper_bound, max_val);
         fail_check = true;
    elseif max_val >= upper_bound % Check if it's still exactly on boundary
        fprintf('DEBUG Violin %d: Nudged max value is STILL >= upper bound (%f): %f\n', i, upper_bound, max_val);
        % This might be okay, but unexpected after nudge
    end
    if any(~isfinite(current_data_clamped))
        fprintf('DEBUG Violin %d: Nudged data contains non-finite values!\n', i);
        fail_check = true;
    end
    if fail_check
         warning('Violin %d: Data nudge/clamping check failed. Potential issues calling ksdensity.', i);
         % error('Violin %d: Nudge/Clamping failed verification.', i);
    end
    % --- End Debug Check ---


    % 4. Call ksdensity
    f_current = []; u_current = []; bb_current = NaN; % Initialize outputs for this iteration
    try
        current_bw_val = []; % Bandwidth for this iteration
        if ~isempty(b)
            current_bw_val = b(i);
            % Check bandwidth validity
            if ~isfinite(current_bw_val) || current_bw_val <= 0
                warning('Violin %d: Invalid bandwidth specified (%f). Using default.', i, current_bw_val);
                current_bw_val = []; % Reset to use default
            end
        end

        if isempty(current_bw_val)
             % Use default bandwidth
             % fprintf('Violin %d: Using default bandwidth.\n', i); % Debug trace
             [f_current, u_current, bb_current]=ksdensity(current_data_clamped, ksdensity_options{:});
        else
             % Use specified valid bandwidth
             % fprintf('Violin %d: Using specified bandwidth: %f\n', i, current_bw_val); % Debug trace
             [f_current, u_current, bb_current]=ksdensity(current_data_clamped, ksdensity_options{:}, 'bandwidth', current_bw_val);
        end

    catch ME
        % --- Catch error to provide more context ---
        fprintf('\n--- ERROR occurred during ksdensity call for Violin %d ---\n', i);
        fprintf('Error message: %s\n', ME.message);
        fprintf('Data characteristics BEFORE ksdensity call:\n');
        fprintf('  Original data size: %d x %d\n', size(Y{i}));
        fprintf('  Number of finite data points: %d\n', numel(current_data));
         if ~isempty(current_data)
            fprintf('  Min/Max of finite data (before clamp): %f / %f\n', min(current_data), max(current_data));
            fprintf('  Min/Max of CLAMPED data passed to ksdensity: %f / %f\n', min(current_data_clamped), max(current_data_clamped));
            fprintf('  Class of clamped data: %s\n', class(current_data_clamped));
            fprintf('  Any non-finite in clamped data? %d\n', any(~isfinite(current_data_clamped)));
            fprintf('  Any clamped < lower_bound (%f)? %d\n', lower_bound, any(current_data_clamped < lower_bound - tolerance));
            fprintf('  Any clamped > upper_bound (%f)? %d\n', upper_bound, any(current_data_clamped > upper_bound + tolerance));
        else
             fprintf('  No finite data points were available.\n');
        end
        fprintf('  ksdensity options used: support = [%f, %f]\n', ksdensity_options{2}(1), ksdensity_options{2}(2));
        if ~isempty(b) && ~isempty(current_bw_val) % Check if specific bw was used
             fprintf('  Specified bandwidth (b(%d)): %f\n', i, current_bw_val);
        else
             fprintf('  Bandwidth: Default (estimated by ksdensity)\n');
        end
        fprintf('  Stack trace:\n');
        for k_err=1:length(ME.stack)
            fprintf('    File: %s, Name: %s, Line: %d\n', ME.stack(k_err).file, ME.stack(k_err).name, ME.stack(k_err).line);
        end
        fprintf('---------------------------------------------------------\n');
        rethrow(ME); % Rethrow the original error after printing info
    end

    % 5. Store bandwidth
    bw(i)=bb_current;

    % 6. Post-process ksdensity results
    f_current(f_current < 0) = 0; % Ensure non-negative density
    max_f = max(f_current);
    if max_f > 0
        f_norm = f_current / max_f * 0.3; % Normalize density for plotting width
    else
        f_norm = zeros(size(f_current)); % Keep it flat if no density
    end

    % 7. Handle potentially different output sizes (u_current)
    if i == 1 || isempty(F) % First non-empty group determines size
        U_rows = numel(u_current);
        [F, U] = deal(NaN(U_rows, nGroups)); % Allocate F and U
        F(:,i) = f_norm;
        U(:,i) = u_current;
    elseif numel(u_current) == U_rows
        % Sizes match, just assign
        F(:,i) = f_norm;
        U(:,i) = u_current;
    else
        % Sizes don't match - INTERPOLATE onto the grid of the first group
        warning('Violin %d: ksdensity returned %d points (expected %d). Interpolating.', i, numel(u_current), U_rows);
        u_common = U(:, find(~all(isnan(U),1), 1, 'first')); % Get u-grid from first valid column
        f_interp = interp1(u_current, f_norm, u_common, 'linear', 0); % Interpolate, fill ends with 0
        F(:,i) = f_interp;
        U(:,i) = u_common; % Use the common grid
    end

end % End of main loop for density calculation

%% Plotting Section
if isempty(F) || all(isnan(F(:)))
    warning('No valid density estimates to plot.');
    return; % Exit if no data could be processed
end

% Determine X coordinates for plotting
if isempty(x)
    x_coords = 1:nGroups;
    setX = 0;
else
    if numel(x) ~= nGroups
        error('Provided ''x'' vector must have the same number of elements as data groups (%d).', nGroups);
    end
    x_coords = x;
    setX = 1;
    if ~isempty(xL)
        disp('Warning: Providing explicit ''x'' coordinates. Ignoring ''xlabel'' strings during axis setup.');
        % Note: Labels can still be set manually after the function call if needed.
    end
end

% Find overall Y range for axis limits, ignoring NaNs
min_U = min(U(:));
max_U = max(U(:));
if isnan(min_U) || isnan(max_U) % Handle case where all data might have been empty/NaN
    min_U = data_support(1) - 0.1*diff(data_support); % Fallback range
    max_U = data_support(2) + 0.1*diff(data_support);
end
y_range = max_U - min_U;
y_lower = min_U - 0.05 * y_range;
y_upper = max_U + 0.05 * y_range;


% Plot the violins
figure_handle = gcf; % Get current figure or create one implicitly
axes_handle = gca; % Get current axes or create one implicitly
hold(axes_handle, 'on');

for i=1:nGroups
    if all(isnan(F(:,i))) || all(isnan(U(:,i)))
        %fprintf('Skipping plotting for group %d due to NaN density/coordinates.\n', i);
        continue; % Skip if no valid density estimate for this group
    end

    % Use the determined x coordinate for this group
    current_x = x_coords(i);

    % Define patch vertices - handle potential NaNs from interpolation/padding
    valid_idx = ~isnan(U(:,i)) & ~isnan(F(:,i));
    u_plot = U(valid_idx, i);
    f_plot = F(valid_idx, i);

    if isempty(u_plot)
        continue; % Skip if no valid points remain
    end

    % Create fill coordinates
    x_fill = [f_plot + current_x; flipud(current_x - f_plot)];
    y_fill = [u_plot; flipud(u_plot)];

    % Assign face color for this group
    current_fc = fc(i,:); % Assumes fc has enough rows now

    % --- Plot Violin Fill ---
    if isempty(lc) || isequal(lower(lc), 'none')
        h(i)=fill(axes_handle, x_fill, y_fill, current_fc, 'FaceAlpha', alp, 'EdgeColor','none');
    else
        h(i)=fill(axes_handle, x_fill, y_fill, current_fc, 'FaceAlpha', alp, 'EdgeColor', lc);
    end
    % Make sure hold is still on
    hold(axes_handle, 'on');

     % --- !!! ADDED: Plot Jittered Data Points !!! ---
    if plotDataPoints && ~isempty(FiniteData{i}) % Check flag and if data exists
        y_points = FiniteData{i};
        n_points = numel(y_points);

        % Calculate jitter amount based on normalized violin half-width (0.3)
        % JitterWidth=1 means potentially span the full 0.3 left/right
        max_jitter = 0.3 * jitterWidth;

        % Generate random offsets for x-coordinates
        % Uniform distribution within [-max_jitter, +max_jitter]
        x_offsets = (rand(n_points, 1) * 2 * max_jitter) - max_jitter;

        % Calculate final x-coordinates for points
        x_data_points = current_x + x_offsets;

        % Plot the points
        scatter(axes_handle, x_data_points, y_points, dataPointSize, ...
                'Marker', dataPointMarker, ...
                'Color', dataPointColor, ...
                'MarkerFaceAlpha', alp); % Use same alpha? Or maybe make configurable?

         % Alternative using plot (simpler styling):
         % plot(axes_handle, x_data_points, y_points, ...
         %      'Marker', dataPointMarker, 'Color', dataPointColor, ...
         %      'MarkerSize', dataPointSize, 'LineStyle', 'none');
    end


    % Plot Mean/Median Lines
    if ~isnan(MX(i)) && plotmean
        % Interpolate width at mean value
        mean_y = MX(i);
        % Need to handle cases where mean is outside the range of u_plot
        if mean_y >= min(u_plot) && mean_y <= max(u_plot)
           width_at_mean_pos = interp1(u_plot, f_plot, mean_y, 'linear', 0); % Width on positive side
           % Handle interpolation returning NaN if mean is exactly at an endpoint without data?
           if isnan(width_at_mean_pos), width_at_mean_pos = 0; end
           x_mean_coords = [current_x - width_at_mean_pos, current_x + width_at_mean_pos];
           p_mean_line = plot(axes_handle, x_mean_coords, [mean_y, mean_y], 'Color', mc, 'LineWidth', 2);
           if i == find(~isnan(MX), 1, 'first'), p(1) = p_mean_line; end % Store handle for legend (only once)
        end
    end
    if ~isnan(MED(i)) && plotmedian
        % Interpolate width at median value
        median_y = MED(i);
        if median_y >= min(u_plot) && median_y <= max(u_plot)
            width_at_median_pos = interp1(u_plot, f_plot, median_y, 'linear', 0); % Width on positive side
            if isnan(width_at_median_pos), width_at_median_pos = 0; end
            x_median_coords = [current_x - width_at_median_pos, current_x + width_at_median_pos];
            p_median_line = plot(axes_handle, x_median_coords, [median_y, median_y], 'Color', medc, 'LineWidth', 2, 'LineStyle', '--'); % Often dashed
            if i == find(~isnan(MED), 1, 'first'), p(2) = p_median_line; end % Store handle for legend (only once)
        end
    end
end % End of plotting loop

hold(axes_handle, 'off');

%% Axis setup and Labels
axis(axes_handle, [min(x_coords)-0.5 max(x_coords)+0.5, y_lower, y_upper]);

set(axes_handle, 'TickLength',[0 0], 'FontSize', 12);
box(axes_handle, 'on');

if setX == 0 % Use default integer ticks if x wasn't specified
    set(axes_handle, 'XTick', x_coords);
    if ~isempty(xL)
        if numel(xL) ~= nGroups
            warning('Number of xlabels does not match number of groups. Using default labels.');
            set(axes_handle, 'XTickLabel', cellstr(num2str(x_coords(:)))); % Default numeric labels
        else
            set(axes_handle, 'XTickLabel', xL);
        end
    else % No labels provided
         set(axes_handle, 'XTickLabel', cellstr(num2str(x_coords(:)))); % Default numeric labels
    end
else % x was specified, use those ticks
     set(axes_handle, 'XTick', x_coords);
% % %      % Labels are usually omitted here or set manually later
% % %      set(axes_handle, 'XTickLabel', []); % Clear default numeric labels if X specified
% % %      % If xL was provided AND setX is 1, we warned earlier. User can set manually.
end


%% Add legend if requested
L = []; % Initialize legend handle
legend_items = gobjects(0);
legend_labels = {};
if plotlegend
    if plotmean && isgraphics(p(1)) % Check if mean handle is valid
        legend_items(end+1) = p(1);
        legend_labels{end+1} = 'Mean';
    end
    if plotmedian && isgraphics(p(2)) % Check if median handle is valid
        legend_items(end+1) = p(2);
        legend_labels{end+1} = 'Median';
    end
    if ~isempty(legend_items)
        L=legend(axes_handle, legend_items, legend_labels);
        set(L,'box','off','FontSize',14, 'Location', 'best') % Changed default location
    end
end

% Return only valid graphics handles for h
h = h(isgraphics(h));

%-------------------------------------------------------------------------
end %of function