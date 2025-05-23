function problem = custom1()
    % This is a toy example to show how to construct a .m file that returns
    % a instance of Problem class representing an unconstrained problem.

    p_struct.name = 'custom1';
    p_struct.x0 = [1; 1];
    p_struct.fun = @fun;

    problem = Problem(p_struct);
end

function f = fun(x)
    f = sum(x.^2);
end