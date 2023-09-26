classdef FeaturedProblem < Problem
    %FEATUREDPROBLEM is optimization problem to be used in the benchmarking with extra features.
    properties (GetAccess = public, SetAccess = private)

        feature
        fun_values
        max_cv_values

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
            obj@Problem(problem.fun, problem.x0, problem.xl, problem.xu, problem.aub, ...
            problem.bub, problem.aeq, problem.beq, problem.cub, problem.ceq);
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
            obj.max_cv_values = [];
        end
        
    end
end