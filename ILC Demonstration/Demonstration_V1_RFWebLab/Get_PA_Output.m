function y = Get_PA_Output(x, coefMat)

[memLen, numCols] = size(coefMat);
memLenM1 = memLen-1;
coefReshaped = reshape(coefMat,1,memLen*numCols);
xLen = length(x);
y = zeros(xLen,1);

degLen = numCols;
for timeIdx = memLen:xLen
    xTerms = zeros(1,memLen*degLen);
    xTime = x(timeIdx-(0:memLenM1));
    xTerms(1,1:memLen) = xTime;
    for degIdx = 2:degLen
        degIdxM1 = degIdx-1;
        startPos = degIdxM1*memLen+1;
        endPos = startPos+memLenM1;
        xTerms(1,startPos:endPos) =                     ...
            xTime.*(abs(xTime).^degIdxM1);
    end
    y(timeIdx) = coefReshaped*xTerms(:);
end
y = transpose(y);
end