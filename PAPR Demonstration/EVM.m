function res = EVM(Err_k, Smax)

N = length(Err_k);
res = (1/Smax) * sqrt((1/N) * sum((abs(Err_k)).^2));

end