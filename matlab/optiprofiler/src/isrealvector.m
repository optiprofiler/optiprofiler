function [isrv, len] = isrealvector(x)
%ISREALVECTOR checks whether x is a real vector. If yes, it returns len = length(x); otherwise, len = NaN.
% N.B.: isrealvector([]) = true
%
%   Function from: https://github.com/zaikunzhang/prima

    if isrealrow(x) || isrealcolumn(x)
        isrv = true;
        len = length(x);
    else
        isrv = false;
        len = NaN;
    end
end