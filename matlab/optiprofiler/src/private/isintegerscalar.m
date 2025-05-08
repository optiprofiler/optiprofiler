function isis = isintegerscalar(x)
%ISINTEGERSCALAR checks whether x is an integer scalar.
% N.B.: isintegerscalar([]) = FALSE, isintegerscalar(NaN) = FALSE, isintegerscalar(inf) = FALSE !!!
%
%   Function from: https://github.com/zaikunzhang/prima

    isis = isrealscalar(x) && (rem(x,1) == 0);

end
