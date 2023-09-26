classdef FeaturedProblem < Problem
    %FEATUREDPROBLEM is optimization problem to be used in the benchmarking with extra features.
    %

    properties (GetAccess = public, SetAccess = private)

        feature
        fun_values
        maxcv_values

    end

    properties (Dependent)

        n_eval

    end

    methods

        function obj = FeaturedProblem(problem, feature)
            %{
            Initialize an optimization problem.

            Parameters
            ----------
            problem : Problem
                Problem to be used in the benchmarking.
            feature : Feature
                Feature to be used in the benchmarking.
            %}
            
            % Copy the problem.
            pb_struct = struct('fun', @(x) problem.fun(x), 'x0', problem.x0, 'xl', problem.xl, ...
                'xu', problem.xu, 'aub', problem.aub, 'bub', problem.bub, 'aeq', problem.aeq, ...
                'beq', problem.beq, 'cub', @(x) problem.cub(x), 'ceq', @(x) problem.ceq(x));
            obj@Problem(pb_struct);
            try
                obj.m_nonlinear_ub = problem.m_nonlinear_ub;
            catch
                % pass
            end
            try
                obj.m_nonlinear_eq = problem.m_nonlinear_eq;
            catch
                % pass
            end

            % Feature to be used in the benchmarking.
            obj.feature = feature;

            % Store the objective function values and maximum constraint violations.
            obj.fun_values = [];
            obj.maxcv_values = [];
        end

        function value = get.n_eval(obj)
            % Return number of objective function evaluations.

            value = length(obj.fun_values);
        end

        function value = get.fun_values(obj)
            value = obj.fun_values;
        end

        function value = get.maxcv_values(obj)
            value = obj.maxcv_values;
        end

        function f = fun(obj, x, seed)
            %{
            Evaluate the objective function.

            Parameters
            ----------
            x : array_like, shape (n,)
                Point at which to evaluate the objective function.
            seed : int, optional
                Seed for the random number generator.

            Returns
            -------
            float
                Value of the objective function at `x`.

            Raises
            ------
            ValueError
                If the argument `x` has an invalid shape.
            %}

            if nargin < 3
                seed = mod(floor(posixtime(datetime('now'))), 2^32);
            end

            % Evaluate the objective function and store the results.
            f_true = fun@Problem(obj, x);
            obj.fun_values = [obj.fun_values, f_true];
            obj.maxcv_values = [obj.maxcv_values, obj.maxcv(x)];

            % Modified the objective function value according to the feature and return the 
            % modified value. We should not store the modified value because the performance 
            % of an optimization solver should be measured using the original objective function.
            f = obj.feature.modifier(x, f_true, seed);
        end
        
    end
end