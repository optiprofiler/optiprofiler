function test_plain(benchmark_id)
    % Test "plain" feature.

    % Go to the directory of this repository.
    cd(fullfile(fileparts(mfilename('fullpath')), '../../../../output'));

    solvers = {@fmincon_test1, @fmincon_test2};
    options.solver_names = {'sqp', 'interior-point'};
    options.feature_name = 'plain';
    options.p_type = 'ubln';
    options.mindim = 11;
    options.maxdim = 11;
    options.benchmark_id = benchmark_id;

    benchmark(solvers, options);
end

function x = fmincon_test1(varargin)
    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'sqp');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end

function x = fmincon_test2(varargin)

    options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'Algorithm', 'interior-point');
    if nargin == 2
        fun = varargin{1};
        x0 = varargin{2};
        x = fmincon(fun, x0, [], [], [], [], [], [], [], options);
    elseif nargin == 4
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        x = fmincon(fun, x0, [], [], [], [], xl, xu, [], options);
    elseif nargin == 8
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, [], options);
    elseif nargin == 10
        fun = varargin{1};
        x0 = varargin{2};
        xl = varargin{3};
        xu = varargin{4};
        aub = varargin{5};
        bub = varargin{6};
        aeq = varargin{7};
        beq = varargin{8};
        cub = varargin{9};
        ceq = varargin{10};
        nonlcon = @(x) deal(cub(x), ceq(x));
        x = fmincon(fun, x0, aub, bub, aeq, beq, xl, xu, nonlcon, options);
    end
end