function problem = custom4()
    % This is a toy example to show how to construct a .m file that returns
    % a instance of Problem class representing a linearly constrained problem.

    p_struct.name = 'custom2';
    p_struct.x0 = [1; 1];
    p_struct.fun = @fun;
    p_struct.xl = [-5; -5];
    p_struct.xu = [5; 5];
    p_struct.aub = [1, 1];
    p_struct.bub = 1;
    p_struct.aeq = [1, -1];
    p_struct.beq = 0;
    p_struct.cub = @cub;
    p_struct.ceq = @ceq;

    problem = Problem(p_struct);
end

function fx = fun(x)
    fx = sum(x.^2);
end

function cubx = cub(x)
    cubx = -sin(x);
end

function ceqx = ceq(x)
    ceqx = sum(x.^4);
end