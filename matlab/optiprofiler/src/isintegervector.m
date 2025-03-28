function isiv = isintegervector(x)
%ISINTEGERVECTOR checks whether x is an integer vector.
% N.B.: isintegervector([]) = FALSE, isintegervector(NaN) = FALSE, isintegervector(inf) = FALSE !!!
%

    isiv = isrealvector(x) && all(rem(x,1) == 0);

end