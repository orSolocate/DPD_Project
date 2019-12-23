function coefMat = Get_Coefficients_Matrix(x, y, memLen, degLen)

xLen = length(x);
memLenM1 = memLen-1;
numCols = degLen;
xTerms = zeros(xLen-memLenM1,memLen*degLen);

for timeIdx = memLen:xLen
    xTime = x(timeIdx-(0:memLenM1));
    xTerms(timeIdx-memLenM1,1:memLen) = xTime;
    for colIdx = 2:degLen
        colIdxM1 = colIdx-1;
        startPos = colIdxM1*memLen+1;
        endPos = startPos+memLen-1;
        xTerms(timeIdx-memLenM1,startPos:endPos) =      ...
            xTime.*abs(xTime).^colIdxM1;
    end
end
xVec = xTerms;
yVec = y(memLen:xLen);
coef = xVec\yVec(:);
coefMat = reshape(coef,memLen,numCols);

end