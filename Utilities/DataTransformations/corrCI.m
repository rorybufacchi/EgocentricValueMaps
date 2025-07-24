function [r_low, r_high] = corrCI(r, n, alpha)
% corrCI calculates the confidence interval for a Pearson correlation coefficient
%
% INPUTS:
%   r     - Pearson correlation coefficient (scalar)
%   n     - Sample size (scalar, must be > 3)
%   alpha - Significance level (e.g., 0.05 for 95% CI)
%
% OUTPUTS:
%   r_low  - Lower bound of the confidence interval
%   r_high - Upper bound of the confidence interval

    if n <= 3
        error('Sample size n must be greater than 3.');
    end
    if abs(r) >= 1
        error('Correlation coefficient must be strictly between -1 and 1.');
    end
    if nargin < 3
        alpha = 0.05; % Default to 95% CI
    end

    % Fisher z-transform
    z = 0.5 * log((1 + r) / (1 - r));
    
    % Standard error of z
    SE = 1 / sqrt(n - 3);
    
    % Critical value from the normal distribution
    z_crit = norminv(1 - alpha/2);
    
    % Confidence interval in z-space
    z_low  = z - z_crit * SE;
    z_high = z + z_crit * SE;
    
    % Inverse Fisher transform
    r_low  = (exp(2*z_low)  - 1) / (exp(2*z_low)  + 1);
    r_high = (exp(2*z_high) - 1) / (exp(2*z_high) + 1);
end
