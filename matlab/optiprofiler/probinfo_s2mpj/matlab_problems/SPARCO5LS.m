function varargout = SPARCO5LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPARCO5LS
%    *********
% 
%    A linear least-squares problem arising in data compression, in the
%    optimization format.
%    The matrix A of the problem is not available explicitly (it is a
%    'spotbox' operator), so only the objective function value (1/2)||r||^2,
%    its gradient J'*r and the product if its Hessian J'*J times v are available.
%    The problem uses the Sparco library tools, and a string containing
%    the path to the master Sparco directory (Sparco-master) MUST be
%    supplied as the last argument in every call to the problem.
%
%    Source:  problem 5 in
%    E. van den Berg, M.~P. Friedlander, G. Hennenfent, F. Herrmann,
%    R. Saab and O. Yilmaz,
%    "Algorithm 890: Sparco: A Testing Framework for Sparse Reconstruction",
%    ACM Transactions on Mathematical Software, vol. 35(19); pp. 1--16, 2009.
%
%    S2MPJ input: Ph. Toint 7 X 2024.
%
%    classification = 'S-CSUR1-AN-2048-0'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;
name = 'SPARCO5LS';
switch(action)
case 'setup'
    nosparco = 1;
    if(nargin>1)
       sparcodir = varargin{end};
       if ( exist( sparcodir ) )
          addpath( sparcodir );
          addpath( [ sparcodir, '/spotbox' ] );
          nosparco = 0;
       end
    end
    if(nosparco)
       disp([ 'ERROR: SPARCO master directory not found.'])
       varargout{1} = NaN;
       if ( nargout > 1 )
          varargout{2} = NaN;
       end
       return
    end
    pbm       = generateProblem(5);
    pbm.name  = name;
    pb.name   = name;
    pb.n      = length( pbm.x0 );
    pb.m      = 0;
    pb.x0     = ones( pb.n, 1 );
%   pb.x0     = pbm.x0;   % solution
    pb.xlower = -Inf*ones(pb.n,1);
    pb.xupper = +Inf*ones(pb.n,1);
    pb.class  = ['S-NLR1-MN-',int2str(pb.n),'-',int2str(pb.m)];
    pbm.objderlvl = 2;
    pb.objderlvl  = 2;
    varargout{1}  = pb;
    varargout{2}  = pbm;
case { 'fx', 'fgx' }     %            ( 'fx', x, sparcodir )
    x  = varargin{1};
    A  = pbm.A;
    b  = pbm.b;
    r  = A*x - b;
    varargout{1} = 0.5 * r' * r;
    if ( nargout > 1 )
       Jt = pbm.A';
       varargout{2} = Jt * r;
    end
case 'fHxv'  %            ( 'fHxv', x, v, sparcodir )
    x  = varargin{1};
    v  = varargin{2};
    J  = pbm.A;
    Jt = pbm.A';
    varargout{1} = Jt*(J*v);
otherwise
    disp([ ' ERROR: action ''',action,''' is unavailable for problem ',name,'.m'])
end

return

end

    