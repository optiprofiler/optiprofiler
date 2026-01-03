function optiprofiler
%OPTIPROFILER is a platform for benchmarking optimization solvers.
%
%   OptiProfiler provides the following functions to facilitate benchmarking
%   of solvers on various test suites.
%
%   Official website:       https://www.optprof.com
%
%   Main function:
%
%   - benchmark: Benchmarks different solvers under different problem libraries
%     with specific features. It generates performance profiles, data profiles,
%     and log-ratio profiles.
%
%   Other tools:
%
%   - s2mpj_load: Loads a problem from the S2MPJ collection.
%
%   - s2mpj_select: Selects problems from the S2MPJ collection based on
%     specified criteria (dimension, constraints, etc.).
%
%   - matcutest_load: Load a problem from the MatCUTEst interface.
%
%   - matcutest_select: Selects problems from the MatCUTEst interface based on
%     specified criteria.
%
%   Try "help benchmark", "help s2mpj_load", etc., for more information.
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Authors:
%               Cunxin HUANG (cun-xin.huang@connect.polyu.hk)
%               Tom M. RAGONNEAU (t.ragonneau@gmail.com)
%               Zaikun ZHANG (zhangzaikun@mail.sysu.edu.cn)
%               Department of Mathematics,
%               Sun Yat-sen University
%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   This file is part of OptiProfiler.