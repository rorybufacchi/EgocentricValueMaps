function CI_r = CalcRhoConfInts(r,N)

% Step 2: Fisher Z-Transformation
Z = 0.5 * log((1 + r) / (1 - r));

% Step 3: Standard Error
SE = 1 / sqrt(N - 3);

% Step 4: Confidence Intervals in Z-space
Z_score = norminv([0.025 0.975]); % For 95% confidence interval
CI_Z = Z + Z_score * SE;

% Step 5: Transform Back to Correlation Space
CI_r = (exp(2*CI_Z) - 1) ./ (exp(2*CI_Z) + 1);

% Display Results
fprintf('Pearson Correlation Coefficient: %f\n', r);
fprintf('95%% Confidence Interval: [%f, %f]\n', CI_r(1), CI_r(2));

end