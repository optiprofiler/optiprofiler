function options = customRuntimeOptions(stage,predecessor)
%CUSTOMRUNTIMEOPTIONS Fresh callback validators bound only to a runtime view.
% The canonical specification retains the original native callbacks. These
% wrappers never enter Feature, saved configuration, or experiment ownership.
    options = stage.options;
    identity = stage.identity;
    n = predecessor.n;
    m_ub = predecessor.m_nonlinear_ub;
    m_eq = predecessor.m_nonlinear_eq;
    if isfield(options,'mod_fun')
        callback = options.mod_fun;
        options.mod_fun = @(x,s,p) scalarOutput(callback,x,s,p,identity);
    end
    if isfield(options,'mod_cub')
        callback = options.mod_cub;
        options.mod_cub = @(x,s,p) vectorOutput(callback(x,s,p),m_ub,identity,'mod_cub');
    end
    if isfield(options,'mod_ceq')
        callback = options.mod_ceq;
        options.mod_ceq = @(x,s,p) vectorOutput(callback(x,s,p),m_eq,identity,'mod_ceq');
    end
    if isfield(options,'mod_x0')
        callback = options.mod_x0;
        options.mod_x0 = @(s,p) vectorOutput(callback(s,p),n,identity,'mod_x0');
    end
    if isfield(options,'mod_bounds')
        callback = options.mod_bounds;
        options.mod_bounds = @(s,p) boundsOutput(callback,s,p,n,identity);
    end
    for key = {'mod_linear_ub','mod_linear_eq'}
        name = key{1};
        if isfield(options,name)
            callback = options.(name);
            options.(name) = @(s,p) linearOutput(callback,s,p,n,identity,name);
        end
    end
    if isfield(options,'mod_affine')
        callback = options.mod_affine;
        options.mod_affine = @(s,p) affineOutput(callback,s,p,n,identity);
    end
end
function value = scalarOutput(callback,x,stream,problem,identity)
    value = callback(x,stream,problem);
    if ~((isnumeric(value) || islogical(value)) && isreal(value) && isscalar(value))
        warning('MATLAB:Feature:InvalidCustomObjective', ...
            'The mod_fun callback of stage %s returned a non-real-scalar value; the observation is NaN.',identity);
        value = NaN;
    else
        value = double(value);
    end
end
function value = vectorOutput(value,count,identity,key)
    if ~((isnumeric(value) || islogical(value)) && isreal(value) ...
            && (isvector(value) || isempty(value)) && numel(value)==count)
        invalid(identity,key,sprintf('a real vector of length %d',count));
    end
    value = double(value(:));
end
function value = matrixOutput(value,shape,identity,key)
    if ~((isnumeric(value) || islogical(value)) && isreal(value) && ismatrix(value) ...
            && isequal(size(value),shape))
        invalid(identity,key,sprintf('a real matrix of shape %d by %d',shape(1),shape(2)));
    end
    value = double(value);
end
function [lower,upper] = boundsOutput(callback,stream,problem,n,identity)
    [lower,upper] = callback(stream,problem);
    lower = vectorOutput(lower,n,identity,'mod_bounds');
    upper = vectorOutput(upper,n,identity,'mod_bounds');
end
function [matrix,rhs] = linearOutput(callback,stream,problem,n,identity,key)
    [matrix,rhs] = callback(stream,problem);
    if isempty(matrix)
        matrix = zeros(0,n);
    else
        matrix = matrixOutput(matrix,[size(matrix,1),n],identity,key);
    end
    rhs = vectorOutput(rhs,size(matrix,1),identity,key);
end
function [matrix,shift,inverse] = affineOutput(callback,stream,problem,n,identity)
    [matrix,shift,inverse] = callback(stream,problem);
    matrix = matrixOutput(matrix,[n,n],identity,'mod_affine');
    shift = vectorOutput(shift,n,identity,'mod_affine');
    inverse = matrixOutput(inverse,[n,n],identity,'mod_affine');
end
function invalid(identity,key,expected)
    error('MATLAB:Feature:InvalidCustomOutput', ...
        'The %s callback of stage %s must return %s.',key,identity,expected);
end
