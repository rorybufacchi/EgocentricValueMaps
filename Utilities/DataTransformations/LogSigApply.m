% Sigmoid Symmetric Transfer Function
function a = LogSigApply(n)
  a = 1 ./ (1 + exp(-n));
end
