function iscs = ischarstr(x)
%ISCHARSTR checks whether an input is a `char` or `string`
%
%   Function from: https://github.com/zaikunzhang/prima

    iscs = (isa(x, 'char') || isa(x, 'string'));
end
