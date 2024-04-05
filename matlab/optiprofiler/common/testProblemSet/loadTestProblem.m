function problem = loadTestProblem(problem_name)
    % This function loads some test problems for plotting the performance profiles (since CUTEst is not available on macOS).

    pb_struct = struct('fun', @boha1, 'x0', ones(2, 1));
    pb1 = Problem(pb_struct);

    pb_struct = struct('fun', @rothyp, 'x0', ones(2, 1));
    pb2 = Problem(pb_struct);

    pb_struct = struct('fun', @rothyp, 'x0', ones(5, 1));
    pb3 = Problem(pb_struct);

    pb_struct = struct('fun', @rothyp, 'x0', ones(10, 1));
    pb4 = Problem(pb_struct);

    pb_struct = struct('fun', @spheref, 'x0', ones(2, 1));
    pb5 = Problem(pb_struct);

    pb_struct = struct('fun', @spheref, 'x0', ones(5, 1));
    pb6 = Problem(pb_struct);

    pb_struct = struct('fun', @spheref, 'x0', ones(10, 1));
    pb7 = Problem(pb_struct);

    pb_struct = struct('fun', @sumpow, 'x0', ones(2, 1));
    pb8 = Problem(pb_struct);

    pb_struct = struct('fun', @sumpow, 'x0', ones(5, 1));
    pb9 = Problem(pb_struct);

    pb_struct = struct('fun', @sumpow, 'x0', ones(10, 1));
    pb10 = Problem(pb_struct);

    pb_struct = struct('fun', @sumsqu, 'x0', ones(2, 1));
    pb11 = Problem(pb_struct);

    pb_struct = struct('fun', @sumsqu, 'x0', ones(5, 1));
    pb12 = Problem(pb_struct);

    pb_struct = struct('fun', @sumsqu, 'x0', ones(10, 1));
    pb13 = Problem(pb_struct);

    pb_struct = struct('fun', @rosen, 'x0', zeros(2, 1));
    pb14 = Problem(pb_struct);

    pb_struct = struct('fun', @rosen, 'x0', zeros(5, 1));
    pb15 = Problem(pb_struct);

    pb_struct = struct('fun', @rosen, 'x0', zeros(10, 1));
    pb16 = Problem(pb_struct);

    switch problem_name
        case 'pb1'
            problem = pb1;
        case 'pb2'
            problem = pb2;
        case 'pb3'
            problem = pb3;
        case 'pb4'
            problem = pb4;
        case 'pb5'
            problem = pb5;
        case 'pb6'
            problem = pb6;
        case 'pb7'
            problem = pb7;
        case 'pb8'
            problem = pb8;
        case 'pb9'
            problem = pb9;
        case 'pb10'
            problem = pb10;
        case 'pb11'
            problem = pb11;
        case 'pb12'
            problem = pb12;
        case 'pb13'
            problem = pb13;
        case 'pb14'
            problem = pb14;
        case 'pb15'
            problem = pb15;
        case 'pb16'
            problem = pb16;
        otherwise
            error('Unknown problem name.');
    end

end
