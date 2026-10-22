function value = referenceTarget(reference, merit_fun, maxcv_init, problem_name, i_run)
%REFERENCETARGET converts one validated raw reference fact into the current merit space.
%   The record was validated when it was constructed. A valid record can still
%   be incompatible with the selected merit: the default merit, for example, is
%   Inf when the violation exceeds max(0.1, 2*maxcv_init). That is reported as
%   an error instead of being replaced by the empirical target, which would
%   silently change the scoring semantics. Used both before solving (preflight)
%   and when scoring, so both report the same message.
    if nargin < 4, problem_name = ''; end
    if nargin < 5, i_run = []; end
    value = merit_fun(reference.fun, reference.maxcv, maxcv_init);
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
        where = '';
        if ~isempty(problem_name)
            where = sprintf(' of problem ''%s''', char(problem_name));
        end
        if ~isempty(i_run)
            where = sprintf('%s (run %d)', where, i_run);
        end
        error('MATLAB:applyReferenceTargets:invalidMerit', ['The supplied reference%s is valid but ' ...
            'incompatible with the selected merit function: the merit function must produce a finite real ' ...
            'scalar for a supplied reference, and it returned %s for kind=''%s'', fun=%.17g, maxcv=%.17g at ' ...
            'maxcv_init=%.17g. The default merit is infinite when the violation exceeds max(0.1, 2 * maxcv_init). ' ...
            'Use a merit_fun that is finite at this reference, or exclude the problem.'], where, ...
            describeValue(value), char(string(reference.kind)), reference.fun, reference.maxcv, maxcv_init);
    end
    value = double(value);
end

function text = describeValue(value)
    if isnumeric(value) && isscalar(value)
        text = num2str(value, 17);
    else
        text = sprintf('a %s of size %s', class(value), mat2str(size(value)));
    end
end
