function f = alwaysError(x)
    % `alwaysError` function is defined for the unit test of `Problem.m`.
    % It activates codes that handles errors when function assignments
    % are made.
   
    error("MyError:alwaysErrorFunction", "This function always throws an error.")
end