function coefMat = System_Identification(yd, u)

% Building the Y matrix, of size NxJ, where J = memLen * degLen and N is
% the legth of yd

memLen = 3;
degLen = 7;

J = memLen * degLen;

Y = zeros(length(yd) - memLen, J);


for n = memLen:length(yd)
    for k = 0:degLen-1
        for m = 0:memLen-1
            curr_col = k * memLen + m + 1;
            Y(n, curr_col) = yd(n-m)*(abs(yd(n-m)))^k;
        end
    end
end

w = inv(ctranspose(Y(3:length(Y), :)) * Y(3:length(Y), :)) * (ctranspose(Y(3:length(Y), :)) * u(3:length(u)));

% Reshaping w to be a matrix instead of a vector

coefMat = reshape(w,memLen,degLen);


end