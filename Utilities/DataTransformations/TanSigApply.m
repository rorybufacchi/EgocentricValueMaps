% Sigmoid Symmetric Transfer Function
function a = TanSigApply(n)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end
