function NanMap(nan_color, main_colormap)
    % nanmap(nan_color, main_colormap) - Function to set a specific color for NaN values.
    %
    % Inputs:
    %   nan_color: a 1x3 RGB triplet that specifies the color for NaN values
    %   main_colormap: the main colormap you are using (e.g., output from jet, hot, etc.)

    % First, set the NaN color in the colormap
    colormap([nan_color; main_colormap]);

end