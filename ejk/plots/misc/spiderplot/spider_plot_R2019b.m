function spider_plot_R2019b(P, options)
%spider_plot_R2019b Create a spider or radar plot with individual axes.
%
% Syntax:
%   spider_plot_R2019b(P)
%   spider_plot_R2019b(P, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   P                - The data points used to plot the spider chart. The
%                      rows are the groups of data and the columns are the
%                      data points. The axes labels and axes limits are
%                      automatically generated if not specified.
%                      [vector | matrix]
%
% Name-Value Pair Arguments:
%   (Optional)
%   AxesLabels       - Used to specify the label each of the axes.
%                      [auto-generated (default) | cell array of character vectors | 'none']
%
%   AxesInterval     - Used to change the number of intervals displayed
%                      between the webs.
%                      [3 (default) | integer]
%
%   AxesPrecision    - Used to change the precision level on the value
%                      displayed on the axes.
%                      [1 (default) | integer | vector]
%
%   AxesDisplay      - Used to change the number of axes in which the
%                      axes text are displayed. 'None' or 'one' can be used
%                      to simplify the plot appearance for normalized data.
%                      ['all' (default) | 'none' | 'one']
%
%   AxesLimits       - Used to manually set the axes limits. A matrix of
%                      2 x size(P, 2). The top row is the minimum axes
%                      limits and the bottow row is the maximum axes limits.
%                      [auto-scaled (default) | matrix]
%
%   FillOption       - Used to toggle color fill option.
%                      ['off' (default) | 'on' | cell array of character vectors]
%
%   FillTransparency - Used to set color fill transparency.
%                      [0.1 (default) | scalar in range (0, 1) | vector]
%
%   Color            - Used to specify the line color, specified as an RGB
%                      triplet. The intensities must be in the range (0, 1).
%                      [MATLAB colors (default) | RGB triplet]
%
%   LineStyle        - Used to change the line style of the plots.
%                      ['-' (default) | '--' | ':' | '-.' | 'none' | cell array of character vectors]
%
%   LineWidth        - Used to change the line width, where 1 point is
%                      1/72 of an inch.
%                      [0.5 (default) | positive value | vector]
%
%   Marker           - Used to change the options.Marker symbol of the plots.
%                      ['o' (default) | 'none' | '*' | 's' | 'd' | ... | cell array of character vectors]
%
%   MarkerSize       - Used to change the options.Marker size, where 1 point is
%                      1/72 of an inch.
%                      [8 (default) | positive value | vector]
%
%   AxesFont         - Used to change the font type of the values
%                      displayed on the axes.
%                      [Helvetica (default) | supported font name]
%
%   LabelFont        - Used to change the font type of the labels.
%                      [Helvetica (default) | supported font name]
%
%   AxesFontSize     - Used to change the font size of the values
%                      displayed on the axes.
%                      [10 (default) | scalar value greater than zero]
%
%   LabelFontSize    - Used to change the font size of the labels.
%                      [10 (default) | scalar value greater than zero]
%
%   Direction        - Used to change the direction of rotation of the
%                      plotted data and axis labels.
%                      ['clockwise' (default) | 'counterclockwise']
%
%   AxesDirection     - Used to change the direction of axes.
%                      ['normal' (default) | 'reverse' | cell array of character vectors]
%
%   AxesLabelsOffset - Used to adjust the position offset of the axes
%                      labels.
%                      [0.1 (default) | positive value]
%
%   AxesScaling      - Used to change the scaling of the axes.
%                      ['linear' (default) | 'log' | cell array of character vectors]
%
%   AxesColor        - Used to change the color of the spider axes.
%                      [grey (default) | RGB triplet | hexadecimal color code]
%
%   AxesLabelsEdge   - Used to change the edge color of the axes labels.
%                      [black (default) | RGB triplet | hexadecimal color code | 'none']
%
%   AxesOffset       - Used to change to axes offset from the origin.
%                      [1 (default) | any integer less than the axes interval]
%
% Examples:
%   % Example 1: Minimal number of arguments. All non-specified, optional
%                arguments are set to their default values. Axes labels
%                and limits are automatically generated and set.
%
%   D1 = [5 3 9 1 2];
%   D2 = [5 8 7 2 9];
%   D3 = [8 2 1 4 6];
%   P = [D1; D2; D3];
%   spider_plot_R2019b(P);
%   legend('D1', 'D2', 'D3', 'Location', 'southoutside');
%
%   % Example 2: Manually setting the axes limits and axes precision.
%                All non-specified, optional arguments are set to their
%                default values.
%
%   D1 = [5 3 9 1 2];
%   D2 = [5 8 7 2 9];
%   D3 = [8 2 1 4 6];
%   P = [D1; D2; D3];
%   spider_plot_R2019b(P,...
%       'AxesLimits', [1, 2, 1, 1, 1; 10, 8, 9, 5, 10],... % [min axes limits; max axes limits]
%       'AxesPrecision', [0, 1, 1, 1, 1]);
%
%   % Example 3: Set fill option on. The fill transparency can be adjusted.
%
%   D1 = [5 3 9 1 2];
%   D2 = [5 8 7 2 9];
%   D3 = [8 2 1 4 6];
%   P = [D1; D2; D3];
%   spider_plot_R2019b(P,...
%       'AxesLabels', {'S1', 'S2', 'S3', 'S4', 'S5'},...
%       'AxesInterval', 2,...
%       'FillOption', {'on', 'on', 'off'},...
%       'FillTransparency', [0.2, 0.1, 0.1]);
%
%   % Example 4: Maximum number of arguments.
%
%   D1 = [5 3 9 1 2];
%   D2 = [5 8 7 2 9];
%   D3 = [8 2 1 4 6];
%   P = [D1; D2; D3];
%   spider_plot_R2019b(P,...
%       'AxesLabels', {'S1', 'S2', 'S3', 'S4', 'S5'},...
%       'AxesInterval', 4,...
%       'AxesPrecision', 0,...
%       'AxesDisplay', 'one',...
%       'AxesLimits', [1, 2, 1, 1, 1; 10, 8, 9, 5, 10],...
%       'FillOption', 'on',...
%       'FillTransparency', 0.2,...
%       'Color', [1, 0, 0; 0, 1, 0; 0, 0, 1],...
%       'LineStyle', {'--', '-', '--'},...
%       'LineWidth', [1, 2, 3],...
%       'Marker', {'o', 'd', 's'},...
%       'MarkerSize', [8, 10, 12],...
%       'AxesFont', 'Times New Roman',...
%       'LabelFont', 'Times New Roman',...
%       'AxesFontSize', 12,...
%       'LabelFontSize', 10,...
%       'Direction', 'clockwise',...
%       'AxesDirection', {'reverse', 'normal', 'normal', 'normal', 'normal'},...
%       'AxesLabelsOffset', 0.1,...
%       'AxesScaling', 'linear',...
%       'AxesColor', [0.6, 0.6, 0.6],...
%       'AxesLabelsEdge', 'none',...
%       'AxesOffset', 1);
%
%   % Example 5: Excel-like radar charts.
%
%   D1 = [5 0 3 4 4];
%   D2 = [2 1 5 5 4];
%   P = [D1; D2];
%   spider_plot_R2019b(P,...
%       'AxesInterval', 5,...
%       'AxesPrecision', 0,...
%       'AxesDisplay', 'one',...
%       'AxesLimits', [0, 0, 0, 0, 0; 5, 5, 5, 5, 5],...
%       'FillOption', 'on',...
%       'FillTransparency', 0.1,...
%       'Color', [139, 0, 0; 240, 128, 128]/255,...
%       'LineWidth', 4,...
%       'Marker', 'none',...
%       'AxesFontSize', 14,...
%       'LabelFontSize', 10,...
%       'AxesColor', [0.8, 0.8, 0.8],...
%       'AxesLabelsEdge', 'none');
%   title('Excel-like Radar Chart',...
%       'FontSize', 14);
%   legend_str = {'D1', 'D2'};
%   legend(legend_str, 'Location', 'southoutside');
%
%   % Example 6: Logarithimic scale on specified axes. Axes limits and axes
%                intervals can be individually set as well.
%
%   D1 = [5 3 9 1 1];
%   D2 = [5 8 7 2 10];
%   D3 = [8 2 1 4 100];
%   P = [D1; D2; D3];
%   spider_plot_R2019b(P,...
%       'AxesInterval', 2,...
%       'AxesPrecision', 0,...
%       'AxesFontSize', 10,...
%       'AxesLabels', {'Linear Scale', 'Linear Scale', 'Linear Scale', 'Linear Scale', 'Logarithimic Scale'},...
%       'AxesScaling', {'linear', 'linear', 'linear', 'linear', 'log'},...
%       'AxesLimits', [1, 1, 1, 1, 1; 10, 10, 10, 10, 100]);
%   legend('D1', 'D2', 'D3', 'Location', 'northeast');
%
%   % Example 7: Spider plot with tiledlayout feature in R2019b.
%
%   D1 = [5 3 9 1 2];
%   D2 = [5 8 7 2 9];
%   D3 = [8 2 1 4 6];
%   P = [D1; D2; D3];
%   t = tiledlayout(2, 2);
%   nexttile;
%   spider_plot_R2019b(P,...
%       'AxesInterval', 1,...
%       'AxesPrecision', 0);
%   nexttile;
%   spider_plot_R2019b(P,...
%       'AxesInterval', 1,...
%       'AxesPrecision', 0);
%   nexttile(3, [1, 2]);
%   spider_plot_R2019b(P,...
%       'AxesInterval', 1,...
%       'AxesPrecision', 0);
%   t.TileSpacing = 'compact';
%   t.Padding = 'compact';
%   title(t, 'Spider Plots');
%
% Author:
%   Moses Yoo, (jyoo at hatci dot com)
%   2020-12-09: Allow fill option and fill transparency for each data group.
%   2020-12-01: Added support for adjust the axes offset from origin.
%   2020-11-30: Allow for one data group without specified axes limits.
%   2020-11-30: Added support for changing axes and label font type.
%   2020-11-06: Fix bug in reverse axes direction feature.
%   2020-10-08: Adjust axes precision to be set to one or more axis.
%   2020-09-30: -Updated examples.
%               -Added feature to change spider axes and axes labels edge color.
%               -Allow logarithmic scale to be set to one or more axis.
%               -Added feature to allow different line styles, line width,
%                marker type, and marker sizes for the data groups.
%               -Allow ability to reverse axes direction.
%   2020-02-12: Fixed condition and added error checking for when only one
%               data group is plotted.
%   2020-01-27: Corrected bug where only 7 entries were allowed in legend.
%   2020-01-06: Added support for tiledlayout feature introduced in R2019b.
%   2019-11-27: Add option to change axes to logarithmic scale.
%   2019-11-15: Add feature to customize the plot rotational direction and
%               the offset position of the axis labels.
%   2019-10-28: Major revision in implementing the new function argument
%               validation feature introduced in R2019b. Replaced previous
%               method of error checking and setting of default values.
%   2019-10-23: Minor revision to set starting axes as the vertical line.
%               Add customization option for font sizes and axes display.
%   2019-10-16: Minor revision to add name-value pairs for customizing
%               color, Marker, and line settings.
%   2019-10-08: Another major revision to convert to name-value pairs and
%               add color fill option.
%   2019-09-17: Major revision to improve speed, clarity, and functionality
%
% Special Thanks:
%   Special thanks to Gabriela Andrade, Andrés Garcia, Jiro Doke,
%   Alex Grenyer, Omar Hadri, Zafar Ali, Christophe Hurlin, Roman,
%   Mariusz Sepczuk, Mohamed Abubakr, Nicolai, Jingwei Too,
%   Cedric Jamet, Richard Ruff & Marie-Kristin Schreiber for their
%   feature recommendations and bug finds.

%%% Argument Validation %%%
arguments
    P (:, :) double
    options.AxesLabels {validateAxesLabels(options.AxesLabels, P)} = cellstr("Label " + (1:size(P, 2)))
    options.AxesInterval (1, 1) double {mustBeInteger, mustBePositive} = 3
    options.AxesPrecision (:, :) double {mustBeInteger, mustBeNonnegative} = 1
    options.AxesDisplay char {mustBeMember(options.AxesDisplay, {'all', 'none', 'one'})} = 'all'
    options.AxesLimits double {validateAxesLimits(options.AxesLimits, P)} = []
    options.FillOption {mustBeMember(options.FillOption, {'off', 'on'})} = 'off'
    options.FillTransparency double {mustBeGreaterThanOrEqual(options.FillTransparency, 0), mustBeLessThanOrEqual(options.FillTransparency, 1)} = 0.1
    options.Color = get(groot,'defaultAxesColorOrder')
    options.LineStyle = '-'
    options.LineWidth (:, :) double {mustBePositive} = 2
    options.Marker = 'o'
    options.MarkerSize (:, :) double {mustBePositive} = 8
    options.AxesFont char = 'Helvetica'
    options.LabelFont char = 'Helvetica'
    options.AxesFontSize (1, 1) double {mustBePositive} = 10
    options.LabelFontSize (1, 1) double {mustBePositive} = 10
    options.Direction char {mustBeMember(options.Direction, {'counterclockwise', 'clockwise'})} = 'clockwise'
    options.AxesDirection = 'normal'
    options.AxesLabelsOffset (1, 1) double {mustBeNonnegative} = 0.1
    options.AxesScaling = 'linear'
    options.AxesColor = [0.6, 0.6, 0.6]
    options.AxesLabelsEdge = 'k'
    options.AxesOffset (1, 1) double {mustBeNonnegative, mustBeInteger} = 1
end

%%% Data Properties %%%
% Point properties
[num_data_groups, num_data_points] = size(P);

%%% Validate Properties %%%
% Check if axes offset is valid
if options.AxesOffset > options.AxesInterval
    error('Error: Invalid axes offset entry. Please enter in an integer value that is between [0, axes_interval].');
end

%%% Validate Axes Precision
% Check if axes precision is numeric
if isnumeric(options.AxesPrecision)
    % Check is length is one
    if length(options.AxesPrecision) == 1
        % Repeat array to number of data points
        options.AxesPrecision = repmat(options.AxesPrecision, num_data_points, 1);
    elseif length(options.AxesPrecision) ~= num_data_points
        error('Error: Please specify the same number of axes precision as number of data points.');
    end
else
    error('Error: Please make sure the axes precision is a numeric value.');
end

%%% Validate Line Style
% Check if line style is a char
if ischar(options.LineStyle)
    % Convert to cell array of char
    options.LineStyle = cellstr(options.LineStyle);
    
    % Repeat cell to number of data groups
    options.LineStyle = repmat(options.LineStyle, num_data_groups, 1);
elseif iscellstr(options.LineStyle) %#ok<*ISCLSTR>
    % Check is length is one
    if length(options.LineStyle) == 1
        % Repeat cell to number of data groups
        options.LineStyle = repmat(options.LineStyle, num_data_groups, 1);
    elseif length(options.LineStyle) ~= num_data_groups
        error('Error: Please specify the same number of line styles as number of data groups.');
    end
else
    error('Error: Please make sure the line style is a char or a cell array of char.');
end

%%% Validate Line Width
% Check if line width is numeric
if isnumeric(options.LineWidth)
    % Check is length is one
    if length(options.LineWidth) == 1
        % Repeat array to number of data groups
        options.LineWidth = repmat(options.LineWidth, num_data_groups, 1);
    elseif length(options.LineWidth) ~= num_data_groups
        error('Error: Please specify the same number of line width as number of data groups.');
    end
else
    error('Error: Please make sure the line width is a numeric value.');
end

%%% Validate Marker
% Check if marker type is a char
if ischar(options.Marker)
    % Convert to cell array of char
    options.Marker = cellstr(options.Marker);
    
    % Repeat cell to number of data groups
    options.Marker = repmat(options.Marker, num_data_groups, 1);
elseif iscellstr(options.Marker)
    % Check is length is one
    if length(options.Marker) == 1
        % Repeat cell to number of data groups
        options.Marker = repmat(options.Marker, num_data_groups, 1);
    elseif length(options.Marker) ~= num_data_groups
        error('Error: Please specify the same number of line styles as number of data groups.');
    end
else
    error('Error: Please make sure the line style is a char or a cell array of char.');
end

%%% Validate Marker Size
% Check if line width is numeric
if isnumeric(options.MarkerSize)
    if length(options.MarkerSize) == 1
        % Repeat array to number of data groups
        options.MarkerSize = repmat(options.MarkerSize, num_data_groups, 1);
    elseif length(options.MarkerSize) ~= num_data_groups
        error('Error: Please specify the same number of line width as number of data groups.');
    end
else
    error('Error: Please make sure the line width is numeric.');
end

%%% Validate Axes Scaling
% Check if axes scaling is valid
if any(~ismember(options.AxesScaling, {'linear', 'log'}))
    error('Error: Invalid axes scaling entry. Please enter in "linear" or "log" to set axes scaling.');
end

% Check if axes scaling is a cell
if iscell(options.AxesScaling)
    % Check is length is one
    if length(options.AxesScaling) == 1
        % Repeat array to number of data points
        options.AxesScaling = repmat(options.AxesScaling, num_data_points, 1);
    elseif length(options.AxesScaling) ~= num_data_points
        error('Error: Please specify the same number of axes scaling as number of data points.');
    end
else
    % Repeat array to number of data points
    options.AxesScaling = repmat({options.AxesScaling}, num_data_points, 1);
end

%%% Validate Axes Direction
% Check if axes direction is a cell
if iscell(options.AxesDirection)
    % Check is length is one
    if length(options.AxesDirection) == 1
        % Repeat array to number of data points
        options.AxesDirection = repmat(options.AxesDirection, num_data_points, 1);
    elseif length(options.AxesDirection) ~= num_data_points
        error('Error: Please specify the same number of axes direction as number of data points.');
    end
else
    % Repeat array to number of data points
    options.AxesDirection = repmat({options.AxesDirection}, num_data_points, 1);
end

%%% Validate Fill Option
% Check if fill option is a cell
if iscell(options.FillOption)
    % Check is length is one
    if length(options.FillOption) == 1
        % Repeat array to number of data groups
        options.FillOption = repmat(options.FillOption, num_data_groups, 1);
    elseif length(options.FillOption) ~= num_data_groups
        error('Error: Please specify the same number of fill option as number of data groups.');
    end
else
    % Repeat array to number of data groups
    options.FillOption = repmat({options.FillOption}, num_data_groups, 1);
end

%%% Validate Fill Transparency
% Check if fill transparency is numeric
if isnumeric(options.FillTransparency)
    % Check is length is one
    if length(options.FillTransparency) == 1
        % Repeat array to number of data groups
        options.FillTransparency = repmat(options.FillTransparency, num_data_groups, 1);
    elseif length(options.FillTransparency) ~= num_data_groups
        error('Error: Please specify the same number of fill transparency as number of data groups.');
    end
else
    error('Error: Please make sure the transparency is a numeric value.');
end

%%% Axes Scaling Properties %%%
% Check axes scaling option
log_index = strcmp(options.AxesScaling, 'log');

% If any log scaling is specified
if any(log_index)
    % Initialize copy
    P_log = P(:, log_index);
    
    % Logarithm of base 10, account for numbers less than 1
    P_log = sign(P_log) .* log10(abs(P_log));
    
    % Minimum and maximun log limits
    min_limit = min(min(fix(P_log)));
    max_limit = max(max(ceil(P_log)));
    recommended_axes_interval = max_limit - min_limit;
    
    % Warning message
    warning('For the log scale values, recommended axes limit is [%i, %i] with an axes interval of %i.',...
        10^min_limit, 10^max_limit, recommended_axes_interval);
    
    % Replace original
    P(:, log_index) = P_log;
end

%%% Figure Properties %%%
% Grab current figure
fig = gcf;

% Set figure background
fig.Color = 'white';

% Reset axes
cla reset;

% Current axes handle
ax = gca;

% Axis limits
hold on;
axis square;
axis([-1, 1, -1, 1] * 1.3);

% Axis properties
ax.XTickLabel = [];
ax.YTickLabel = [];
ax.XColor = 'none';
ax.YColor = 'none';

% Polar increments
theta_increment = 2*pi/num_data_points;
full_interval = options.AxesInterval+1;
rho_offset = options.AxesOffset/full_interval;

%%% Scale Data %%%
% Pre-allocation
P_scaled = zeros(size(P));
axes_range = zeros(3, num_data_points);

% Check axes scaling option
axes_direction_index = strcmp(options.AxesDirection, 'reverse');

% Iterate through number of data points
for ii = 1:num_data_points
    % Check for one data group and no axes limits
    if num_data_groups == 1 && isempty(options.AxesLimits)
        % Group of points
        group_points = P(:, :);
    else
        % Group of points
        group_points = P(:, ii);
    end
    
    % Check for log axes scaling option
    if log_index(ii)
        % Minimum and maximun log limits
        min_value = min(fix(group_points));
        max_value = max(ceil(group_points));
    else
        % Automatically the range of each group
        min_value = min(group_points);
        max_value = max(group_points);
    end
    
    % Range of min and max values
    range = max_value - min_value;
    
    % Check if options.AxesLimits is not empty
    if ~isempty(options.AxesLimits)
        % Check for log axes scaling option
        if log_index(ii)
            % Logarithm of base 10, account for numbers less than 1
            options.AxesLimits(:, ii) = sign(options.AxesLimits(:, ii)) .* log10(abs(options.AxesLimits(:, ii)));
        end
            
        % Manually set the range of each group
        min_value = options.AxesLimits(1, ii);
        max_value = options.AxesLimits(2, ii);
        range = max_value - min_value;
        
        % Check if the axes limits are within range of points
        if min_value > min(group_points) || max_value < max(group_points)
            error('Error: Please make the manually specified axes limits are within range of the data points.');
        end
    end
    
    % Scale points to range from [0, 1]
    P_scaled(:, ii) = ((P(:, ii) - min_value) / range);
    
    % If reverse axes direction is specified
    if axes_direction_index(ii)
        % Store to array
        axes_range(:, ii) = [max_value; min_value; range];
        P_scaled(:, ii) = -(P_scaled(:, ii) - 1);
    else
        % Store to array
        axes_range(:, ii) = [min_value; max_value; range];
    end
    
    % Add offset of [rho_offset] and scaling factor of [1 - rho_offset]
    P_scaled(:, ii) = P_scaled(:, ii) * (1 - rho_offset) + rho_offset;
end

%%% Polar Axes %%%
% Polar coordinates
rho_increment = 1/full_interval;
rho = 0:rho_increment:1;

% Check specified direction of rotation
switch options.Direction
    case 'counterclockwise'
        % Shift by pi/2 to set starting axis the vertical line
        theta = (0:theta_increment:2*pi) + (pi/2);
    case 'clockwise'
        % Shift by pi/2 to set starting axis the vertical line
        theta = (0:-theta_increment:-2*pi) + (pi/2);
end

% Remainder after using a modulus of 2*pi
theta = mod(theta, 2*pi);

% Iterate through each theta
for ii = 1:length(theta)-1
    % Convert polar to cartesian coordinates
    [x_axes, y_axes] = pol2cart(theta(ii), rho);
    
    % Plot webs
    h = plot(x_axes, y_axes,...
        'LineWidth', 1.5,...
        'Color', options.AxesColor);
    
    % Turn off legend annotation
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% Iterate through each rho
for ii = 2:length(rho)
    % Convert polar to cartesian coordinates
    [x_axes, y_axes] = pol2cart(theta, rho(ii));
    
    % Plot axes
    h = plot(x_axes, y_axes,...
        'Color', options.AxesColor);
    
    % Turn off legend annotation
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

% Set end index depending on axes display argument
switch options.AxesDisplay
    case 'all'
        theta_end_index = length(theta)-1;
    case 'one'
        theta_end_index = 1;
    case 'none'
        theta_end_index = 0;
end

% Rho start index and offset interval
rho_start_index = options.AxesOffset+1;
offset_interval = full_interval - options.AxesOffset;

% Iterate through each theta
for ii = 1:theta_end_index
    % Convert polar to cartesian coordinates
    [x_axes, y_axes] = pol2cart(theta(ii), rho);
    
    % Iterate through points on isocurve
    for jj = rho_start_index:length(rho)
        % Axes increment range
        min_value = axes_range(1, ii);
        range = axes_range(3, ii);
        
        % If reverse axes direction is specified
        if axes_direction_index(ii)
            % Axes increment value
            axes_value = min_value - (range/offset_interval) * (jj-rho_start_index);
        else
            % Axes increment value
            axes_value = min_value + (range/offset_interval) * (jj-rho_start_index);
        end
        
        % Check for log axes scaling option
        if log_index(ii)
            % Exponent to the tenth power
            axes_value = 10^axes_value;
        end
        
        % Display axes text
        text_str = sprintf(sprintf('%%.%if', options.AxesPrecision(ii)), axes_value);
        text(x_axes(jj), y_axes(jj), text_str,...
            'Units', 'Data',...
            'Color', 'k',...
            'FontName', options.AxesFont,...
            'FontSize', options.AxesFontSize,...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle');
    end
end

%%% Plot %%%
% Fill option index
fill_option_index = strcmp(options.FillOption, 'on');

% Iterate through number of data groups
for ii = 1:num_data_groups
    % Convert polar to cartesian coordinates
    [x_points, y_points] = pol2cart(theta(1:end-1), P_scaled(ii, :));
    
    % Make points circular
    x_circular = [x_points, x_points(1)];
    y_circular = [y_points, y_points(1)];
    
    % Plot data points
    plot(x_circular, y_circular,...
        'LineStyle', options.LineStyle{ii},...
        'Marker', options.Marker{ii},...
        'Color', options.Color(ii, :),...
        'LineWidth', options.LineWidth(ii),...
        'MarkerSize', options.MarkerSize(ii),...
        'MarkerFaceColor', options.Color(ii, :));
    
    % Check if fill option is toggled on
    if fill_option_index(ii)
        % Fill area within polygon
        h = patch(x_circular, y_circular, options.Color(ii, :),...
            'EdgeColor', 'none',...
            'FaceAlpha', options.FillTransparency(ii));
        
        % Turn off legend annotation
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

% Get object type of children
children_type = get(fig.Children, 'Type');

% Check if children is titled layout
if strcmp(children_type, 'titledlayout')
    % Set handle to axes children
    obj = ax.Children;
else
    % Set handle to axes
    obj = ax;
end

% Find object handles
text_handles = findobj(obj.Children,...
    'Type', 'Text');
patch_handles = findobj(obj.Children,...
    'Type', 'Patch');
isocurve_handles = findobj(obj.Children,...
    'Color', options.AxesColor,...
    '-and', 'Type', 'Line');
plot_handles = findobj(obj.Children, '-not',...
    'Color', options.AxesColor,...
    '-and', 'Type', 'Line');

% Manually set the stack order
uistack(text_handles, 'bottom');
uistack(plot_handles, 'bottom');
uistack(patch_handles, 'bottom');
uistack(isocurve_handles, 'bottom');

%%% Labels %%%
% Check labels argument
if ~strcmp(options.AxesLabels, 'none')
    % Convert polar to cartesian coordinates
    [x_axes, y_axes] = pol2cart(theta, rho(end));
    
    % Iterate through number of data points
    for ii = 1:length(options.AxesLabels)
        % Angle of point in radians
        theta_point = theta(ii);
        
        % Find out which quadrant the point is in
        if theta_point == 0
            quadrant = 0;
        elseif theta_point == pi/2
            quadrant = 1.5;
        elseif theta_point == pi
            quadrant = 2.5;
        elseif theta_point == 3*pi/2
            quadrant = 3.5;
        elseif theta_point == 2*pi
            quadrant = 0;
        elseif theta_point > 0 && theta_point < pi/2
            quadrant = 1;
        elseif theta_point > pi/2 && theta_point < pi
            quadrant = 2;
        elseif theta_point > pi && theta_point < 3*pi/2
            quadrant = 3;
        elseif theta_point > 3*pi/2 && theta_point < 2*pi
            quadrant = 4;
        end
        
        % Adjust label alignment depending on quadrant
        switch quadrant
            case 0
                horz_align = 'left';
                vert_align = 'middle';
                x_pos = options.AxesLabelsOffset;
                y_pos = 0;
            case 1
                horz_align = 'left';
                vert_align = 'bottom';
                x_pos = options.AxesLabelsOffset;
                y_pos = options.AxesLabelsOffset;
            case 1.5
                horz_align = 'center';
                vert_align = 'bottom';
                x_pos = 0;
                y_pos = options.AxesLabelsOffset;
            case 2
                horz_align = 'right';
                vert_align = 'bottom';
                x_pos = -options.AxesLabelsOffset;
                y_pos = options.AxesLabelsOffset;
            case 2.5
                horz_align = 'right';
                vert_align = 'middle';
                x_pos = -options.AxesLabelsOffset;
                y_pos = 0;
            case 3
                horz_align = 'right';
                vert_align = 'top';
                x_pos = -options.AxesLabelsOffset;
                y_pos = -options.AxesLabelsOffset;
            case 3.5
                horz_align = 'center';
                vert_align = 'top';
                x_pos = 0;
                y_pos = -options.AxesLabelsOffset;
            case 4
                horz_align = 'left';
                vert_align = 'top';
                x_pos = options.AxesLabelsOffset;
                y_pos = -options.AxesLabelsOffset;
        end
        
        % Display text label
        text(x_axes(ii)+x_pos, y_axes(ii)+y_pos, options.AxesLabels{ii},...
            'Units', 'Data',...
            'HorizontalAlignment', horz_align,...
            'VerticalAlignment', vert_align,...
            'EdgeColor', options.AxesLabelsEdge,...
            'BackgroundColor', 'w',...
            'FontName', options.LabelFont,...
            'FontSize', options.LabelFontSize);
    end
end
end

%%% Custom Validation Functions %%%
% Validate axes labels
function validateAxesLabels(axes_labels, P)
if ~isequal(axes_labels, 'none')
    validateattributes(axes_labels, {'cell'}, {'size', [1, size(P, 2)]}, mfilename, 'axes_labels')
end
end

% Validate axes limits
function validateAxesLimits(axes_limits, P)
if ~isempty(axes_limits)
    validateattributes(axes_limits, {'double'}, {'size', [2, size(P, 2)]}, mfilename, 'axes_limits')
    
    % Lower and upper limits
    lower_limits = axes_limits(1, :);
    upper_limits = axes_limits(2, :);
    
    % Difference in upper and lower limits
    diff_limits = upper_limits - lower_limits;
    
    % Check to make sure upper limit is greater than lower limit
    if any(diff_limits < 0)
        error('Error: Please make sure max axes limits are greater than the min axes limits.');
    end
    
    % Check the range of axes limits
    if any(diff_limits == 0)
        error('Error: Please make sure the min and max axes limits are different.');
    end
end
end