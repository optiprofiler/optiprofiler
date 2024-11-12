function varargout = SPARCO2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPARCO2
%    *********
% 
%    A linear least-squares problem arising in data compression.
%    The matrix A of the problem is not available explicitly (it is a
%    'spotbox' operator), so only the residual values ('cx') and the products
%    A times v ('cJxv) and A' times  v ('cJtxv') are available.
%    The problem uses the Sparco library tools, and a string containing
%    the path to the master Sparco directory (Sparco-master) MUST be
%    supplied as the last argument in every call to the problem.
%
%    Source:  problem 2 in
%    E. van den Berg, M.~P. Friedlander, G. Hennenfent, F. Herrmann,
%    R. Saab and O. Yilmaz,
%    "Algorithm 890: Sparco: A Testing Framework for Sparse Reconstruction",
%    ACM Transactions on Mathematical Software, vol. 35(19); pp. 1--16, 2009.
%
%    S2MPJ input: Ph. Toint 4 X 2024.
%
%    classification = 'S-CNLR1-AN-1024-1024'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;
name = 'SPARCO2';
switch(action)
case 'setup'
    sparcodir = varargin{end};
    if ( exist( sparcodir ) )
       addpath( sparcodir );
       addpath( [ sparcodir, '/spotbox' ] );
    else
       disp([ 'ERROR: SPARCO master directory ',sparcodir,' does not exist.'])
       varargout{1} = NaN;
       return
    end
    pbm       = generateProblem(2);
    pbm.name  = name;
    pb.name   = name;
    pb.n      = length( pbm.x0 );
    pb.m      = length( pbm.b );
    pb.x0     = ones( pb.n, 1);
    pb.xlower = -Inf*ones(pb.n,1);
    pb.xupper = +Inf*ones(pb.n,1);
    pb.class  = ['S-NLR1-MN-',int2str(pb.n),'-',int2str(pb.m)];
    pbm.conderlvl = 2;
    pb.conderlvl  = 2;
    varargout{1}  = pb;
    varargout{2}  = pbm;
case 'cx'     %            ( 'cx', x, sparcodir )
    x  = varargin{1};
    A  = pbm.A;
    b  = pbm.b;
    varargout{1} = A*x - b;
case 'cJxv'   %            ( 'cJxv', x, v, sparcodir )
    x  = varargin{1};
    v  = varargin{2};
    x  = ones(length(x),1);
    J  = pbm.A;
    varargout{1} = J*v;
case 'cJtxv'  %            ( 'cJtxv', x, v, sparcodir )
    x  = varargin{1};
    v  = varargin{2};
    Jt = pbm.A';
    varargout{1} = Jt*v;
otherwise
    disp([ ' ERROR: action ',action,' is unavailable for problem ',name,'.m'])
    varargout{1} = NaN;
end

return

end

