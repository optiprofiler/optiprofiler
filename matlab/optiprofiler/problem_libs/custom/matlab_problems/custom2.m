function problem = custom2()
    % This is a toy example to show how to construct a .m file that returns
    % a instance of Problem class representing a bound-constrained problem.

    p_struct.name = 'custom2';
    p_struct.x0 = [1; 1];
    p_struct.fun = @fun;
    p_struct.xl = [-5; -5];
    p_struct.xu = [5; 5];

    problem = Problem(p_struct);
end

function f = fun(x)
    f = sum(x.^2);
end