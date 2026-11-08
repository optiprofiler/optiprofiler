function [entry, rule] = invalidConstraintEntry(value, name)
%INVALIDCONSTRAINTENTRY First entry of VALUE without a meaning as the constraint data NAME.
% ENTRY is empty if there is none, and RULE states what NAME must hold. Such
% entries used to be measured as if they had a meaning: xl = NaN counted as a
% bound that holds at every point, xl = Inf as no bound, bub = -Inf as no
% constraint, and a coefficient of -Inf satisfies its row at every positive
% point. Violations of finite data, such as xl > xu, are measured and are not
% refused here.
    entries = value(:);
    if issparse(value)
        entries = nonzeros(value);
    end
    switch name
        case 'xl'
            valid = isfinite(entries) | entries == -Inf;
            rule = 'a lower bound must be finite, or -Inf for no bound';
        case 'xu'
            valid = isfinite(entries) | entries == Inf;
            rule = 'an upper bound must be finite, or Inf for no bound';
        case 'bub'
            valid = isfinite(entries) | entries == Inf;
            rule = 'a right-hand side must be finite, or Inf for no constraint';
        case 'beq'
            valid = isfinite(entries);
            rule = 'an equality right-hand side must be finite';
        otherwise
            valid = isfinite(entries);
            rule = 'a coefficient must be finite';
    end
    entry = entries(find(~valid, 1));
end
