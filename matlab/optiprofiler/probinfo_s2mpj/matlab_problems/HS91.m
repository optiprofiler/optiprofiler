function varargout = HS91(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS91
%    *********
% 
%    A time-optimal heat conduction problem.
% 
%    Source: problem 91 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, September 1991.
%      SAVEs removed December 3rd 2014
% 
%    classification = 'C-CQOR2-MN-5-1'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS91';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 5;
        v_('EPS') = 0.01;
        v_('EPSSQR') = v_('EPS')*v_('EPS');
        v_('-EPSSQR') = -1.0*v_('EPSSQR');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CON',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('CON')) = v_('-EPSSQR');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('2'):v_('N')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = 0.5;
            end
        end
        for I=v_('2'):v_('2'):v_('N')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = -0.5;
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = -0.5;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eH',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['O',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'H';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eH';
        ielftype(ie) = iet_('eH');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['O',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_('CON');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('H');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL HS91SOL              1.35919D+00
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MN-5-1';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0;
                varargout{3} = H_;
            end
        end

    case 'eH'

        EV_  = varargin{1};
        iel_ = varargin{2};
        [F,g_,H_] = eval91(EV_);
        varargout{1} = F;
        if(nargout>1)
            varargout{2} = g_;
            if(nargout>2)
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g, H ] = eval91( x )

%   A translation of the Fortran code present in the SIF file.

    n     = length(x);
    g     = zeros(n,1);
    H     = sparse(n,n);
    A     = zeros(30,1);
    R     = zeros(30,30);
    S     = zeros(30,1);
    RHO   = zeros(30,1);
    DRHO  = zeros(30,n);
    D2RHO = zeros(30,n,n);
    P     = zeros(n+1,1);
    
    mu = [ 8.6033358901938017e-01,  3.4256184594817283e+00, 6.4372981791719468e+00,  9.5293344053619631e+00, ...
           1.2645287223856643e+01,  1.5771284874815882e+01, 1.8902409956860023e+01,  2.2036496727938566e+01, ...
           2.5172446326646664e+01,  2.8309642854452012e+01, 3.1447714637546234e+01,  3.4586424215288922e+01, ...
           3.7725612827776501e+01,  4.0865170330488070e+01, 4.4005017920830845e+01,  4.7145097736761031e+01, ...
           5.0285366337773652e+01,  5.3425790477394663e+01, 5.6566344279821521e+01,  5.9707007305335459e+01, ...
           6.2847763194454451e+01,  6.5988598698490392e+01, 6.9129502973895256e+01,  7.2270467060308960e+01, ...
           7.5411483488848148e+01,  7.8552545984242926e+01, 8.1693649235601683e+01,  8.4834788718042290e+01, ...
           8.7975960552493220e+01,  9.1117161394464745e+01 ];

    T = 2.0 / 15.0;
    for i = 1:30
        MUI    = mu(i);
        SMUI   = sin(MUI);
        CMUI   = cos(MUI);
        AI     = 2.0*SMUI/(MUI+SMUI*CMUI);
        A( i ) = AI;
        S( i ) = 2.0*AI*(CMUI-SMUI/MUI);
        AIMUI2 = AI * MUI^2;
        for j = 1:i
            if ( i ~= j ) 
                MUJ    = mu(j);
                R(i,j) = 0.5*(sin(MUI+MUJ )/(MUI+MUJ)+sin(MUI-MUJ )/(MUI-MUJ))*AIMUI2*A(j)*MUJ^2;
                R(j,i) = R(i,j);
            else
                R(i,i) = 0.5*(1.0+0.5*sin(MUI+MUI)/MUI)*AIMUI2^2;
            end
        end
    end

%                                  n   2
%  Calculate the functions p(x) = SUM x .
%                           j     i=j  i

    for k = n:-1:1
        P(k) = P(k+1)+x(k)^2;
    end

%  Calculate the functions rho.

    for j = 1:30
        MUJ2 = mu(j)*mu(j);
        U    = exp(-MUJ2*P(1));
        for k =1:n
            DRHO(j,k) = 2.0*U*x(k);
            for l = k:n
                D2RHO(j,k,l) = -4.0*MUJ2*U*x(k)*x(l);
                if ( l == k )
                    D2RHO(j,k,l) = D2RHO(j,k,l)+2.0*U;
                end
            end
        end
        ALPHA = -2.0;
        for i = 2:n
            EU = ALPHA*exp(-MUJ2*P(i));
            U  = U+EU;
            for k = i:n
                DRHO(j,k) = DRHO(j,k)+2.0*EU*x(k);
                for l = k:n
                    D2RHO(j,k,l) = D2RHO(j,k,l)-4.0*MUJ2*EU*x(k)*x(l);
                    if ( l == k )
                         D2RHO(j,k,l) = D2RHO(j,k,l)+2.0*EU;
                    end
                end
            end
            ALPHA = - ALPHA;
        end
        U      = U+0.5*ALPHA;
        RHO(j) = -U/MUJ2;
    end

%  Evaluate the function and derivatives.

    f = T;
    for i = 1:30
        SI   = S(i);
        RHOI = RHO(i);
        f    = f+SI*RHOI;
        for k = 1:n
            g(k) = g(k)+SI*DRHO(i,k);
            for l = k:n
               H(k,l) = H(k,l)+SI*D2RHO(i,k,l);
            end
        end
        for j = 1:30
            RIJ  = R(i,j);
            RHOJ = RHO(j);
            f    = f+RIJ*RHOI*RHOJ;
            for k = 1:n
                g(k) = g(k)+RIJ*(RHOI*DRHO(j,k)+RHOJ*DRHO(i,k));
                for l = k:n
                    H(k,l) = H(k,l)+RIJ*(RHOI*D2RHO(j,k,l)+RHOJ*D2RHO(i,k,l)+DRHO(i,k)*DRHO(j,l)+DRHO(j,k)*DRHO(i,l));
                end
            end
        end
    end

%   Symmetrize the Hessian.

    for k = 1:n
        for l = k+1:n
            H(l,k) = H(k,l);
        end
    end

return

end
