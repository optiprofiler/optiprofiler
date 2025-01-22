function isrs = isrealscalar(x)
%ISREALSCALAR checks whether x is a real scalar.
% N.B.: isrealscalar([]) = FALSE, isrealscalar(NaN) = TRUE, isrealscalar(inf) = TRUE!!!
%
%   Function from: https://github.com/zaikunzhang/prima

    isrs = isnumeric(x) && isreal(x) && isscalar(x);
end
