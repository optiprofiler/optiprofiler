function problem = toy_example()
    % This is a toy example to show how to construct a .m file that returns
    % a instance of Problem class.
    
    p_struct.x0 = [1; 1];
    p_struct.fun = @obj;

    problem = Problem(p_struct);
end

function f = obj(x)
    f = sum(x.^2);
end