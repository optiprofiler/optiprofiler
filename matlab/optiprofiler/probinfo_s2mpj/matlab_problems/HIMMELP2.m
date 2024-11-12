function varargout = HIMMELP2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELP2
%    *********
% 
%    A nonlinear problem with inequality constraints, attributed to Himmelblau
%    by B.N. Pshenichnyj (case I)
% 
%    The problem is nonconvex and has redundant constraints at the solution.
% 
%    Source: 
%    B.N. Pshenichnyj
%    "The Linearization Method for Constrained Optimization",
%    Springer Verlag, SCM Series 22, Heidelberg, 1994
% 
%    SIF input: Ph. Toint, December 1994.
% 
%    classification = 'C-COQR2-AN-2-1'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELP2';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
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
        v_('B1') = 0.1963666677;
        v_('B1') = 75.0+v_('B1');
        v_('B2') = -.8112755343;
        v_('B2') = -3.0+v_('B2');
        v_('B6') = -.8306567613;
        v_('B6') = -6.0+v_('B6');
        v_('-B2') = -1.0*v_('B2');
        v_('-B6') = -1.0*v_('B6');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-B2')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-B2');
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('-B6')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('-B6');
        end
        [ig,ig_] = s2mpjlib('ii','C',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C';
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
        pbm.gconst(ig_('OBJ')) = v_('B1');
        pbm.gconst(ig_('C')) = 700.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.0;
        pb.xupper(ix_('X1')) = 95.0;
        pb.xlower(ix_('X2'),1) = 0.0;
        pb.xupper(ix_('X2')) = 75.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('X1'),1) = 95.0;
        pb.x0(ix_('X2'),1) = 10.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOBNL',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'OB';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOBNL';
        ielftype(ie) = iet_('eOBNL');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X1X2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OB');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X1X2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                -62.053869846
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COQR2-AN-2-1';
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

    case 'eOBNL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B3 = .1269366345;
        B4 = -0.20567665;
        B4 = 0.01*B4;
        B5 = 0.103450e-4;
        B7 = .0302344793;
        B8 = -0.12813448;
        B8 = 0.01*B8;
        B9 = 0.352599e-4;
        B10 = -0.2266e-6;
        B11 = 0.2564581253;
        B12 = -.003460403;
        B13 = 0.135139e-4;
        B14 = -.1064434908;
        B14 = B14-28.0;
        B15 = -0.52375e-5;
        B16 = -0.63e-8;
        B17 = 0.7e-9;
        B18 = 0.3405462;
        B18 = 0.001*B18;
        B19 = -0.16638e-5;
        B20 = -2.86731123;
        B20 = B20-0.92e-8;
        A = B7*EV_(1)+B8*EV_(1)^2+B9*EV_(1)^3+B10*EV_(1)^4;
        DADX = B7+2.0*B8*EV_(1)+3.0*B9*EV_(1)^2+4.0*B10*EV_(1)^3;
        D2ADXX = 2.0*B8+6.0*B9*EV_(1)+12.0*B10*EV_(1)^2;
        B = B18*EV_(1)+B15*EV_(1)^2+B16*EV_(1)^3;
        DBDX = B18+2.0*B15*EV_(1)+3.0*B16*EV_(1)^2;
        D2BDXX = 2.0*B15+6.0*B16*EV_(1);
        C = B3*EV_(1)^2+B4*EV_(1)^3+B5*EV_(1)^4;
        DCDX = 2.0*B3*EV_(1)+3.0*B4*EV_(1)^2+4.0*B5*EV_(1)^3;
        D2CDXX = 2.0*B3+6.0*B4*EV_(1)+12.0*B5*EV_(1)^2;
        F = B11*EV_(2)^2+B12*EV_(2)^3+B13*EV_(2)^4;
        DFDY = 2.0*B11*EV_(2)+3.0*B12*EV_(2)^2+4.0*B13*EV_(2)^3;
        D2FDYY = 2.0*B11+6.0*B12*EV_(2)+12.0*B13*EV_(2)^2;
        G = B17*EV_(1)^3+B19*EV_(1);
        DGDX = B19+3.0*B17*EV_(1)^2;
        D2GDXX = 6.0*B17*EV_(1);
        E = exp(0.0005*EV_(1)*EV_(2));
        DEDX = 0.0005*EV_(2)*E;
        DEDY = 0.0005*EV_(1)*E;
        D2EDXX = 0.0005*EV_(2)*DEDX;
        D2EDXY = 0.0005*(EV_(2)*DEDY+E);
        D2EDYY = 0.0005*EV_(1)*DEDY;
        varargout{1} = C+EV_(2)*A+F+B14/(1.0+EV_(2))+B*EV_(2)^2+G*EV_(2)^3+B20*E;
        if(nargout>1)
            g_(1,1) = DCDX+EV_(2)*DADX+DBDX*EV_(2)^2+DGDX*EV_(2)^3+B20*DEDX;
            g_(2,1) = A+DFDY-B14/(1.0+EV_(2))^2+2.0*B*EV_(2)+3.0*G*EV_(2)^2+B20*DEDY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = D2CDXX+EV_(2)*D2ADXX+D2BDXX*EV_(2)^2+D2GDXX*EV_(2)^3+B20*D2EDXX;
                H_(1,2) = DADX+2.0*EV_(2)*DBDX+3.0*DGDX*EV_(2)^2+B20*D2EDXY;
                H_(2,1) = H_(1,2);
                H_(2,2) = D2FDYY+2.0*B14/(1.0+EV_(2))^3+2.0*B+6.0*G*EV_(2)+B20*D2EDYY;
                varargout{3} = H_;
            end
        end

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
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

