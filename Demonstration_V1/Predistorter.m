function [u_star, buffer] = Predistorter(u, buffer, coefMat, g)

g_avg = g;
yd = u * g_avg;

[memLen, degLen] = size(coefMat);

if (length(buffer)<(memLen-1))
    u_star = u;
    buffer = [yd buffer];
else
    yMat = zeros(memLen, degLen);
    for k=1:degLen
        yMat(1,k) = yd * abs(yd)^(k-1);
        for m=1:(memLen-1)
            yMat(m+1,k) = buffer(m) * abs(buffer(m))^(k-1);
        end
    end
    u_star = sum(sum(yMat .* coefMat));
    buffer = [yd buffer(1:memLen-1)];
end
end