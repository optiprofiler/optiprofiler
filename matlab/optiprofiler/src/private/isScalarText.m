function valid = isScalarText(value)
%ISSCALARTEXT One piece of text: a char row or a string scalar, not a string array.
    valid = (ischar(value) && (isrow(value) || isempty(value))) || ...
        (isstring(value) && isscalar(value) && ~ismissing(value));
end
