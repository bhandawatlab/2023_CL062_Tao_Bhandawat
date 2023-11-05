function [out] = round2NearestInterval(A,interval,numPow,roundup)
if isempty(interval)
    interval = 5;
end
if isempty(numPow)
    numPow = floor(log10(abs(A)));
end
if roundup
    out = ceil(A./((10.^(numPow-1)).*interval)).*((10.^(numPow-1)).*interval);
else
    out = floor(A./((10.^(numPow-1)).*interval)).*((10.^(numPow-1)).*interval);
end
end