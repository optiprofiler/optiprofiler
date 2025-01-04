function varargout = LUKVLI6(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKVLI6
%    *********
% 
%    Source: Problem 5.6, Generalized Broyden banded function with 
%    exponential constraints, due to L. Luksan and J. Vlcek,
%    "Sparse and partially separable test problems for 
%    unconstrained and equality constrained optimization",
%    Technical Report 767, Inst. Computer Science, Academy of Sciences
%    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
% 
%    Equality constraints changed to inequalities
% 
%    SIF input: Nick Gould, April 2001
% 
%    classification = 'C-COOR2-AY-V-V'
% 
%    some useful parameters, including N, the number of variables.
% 
%       Alternative values for the SIF file parameters:
% IE N                   99             $-PARAMETER
% IE N                   999            $-PARAMETER
% IE N                   9999           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKVLI6';

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
        if(nargs<1)
            v_('N') = 9;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   99999          $-PARAMETER
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N+1') = 1+v_('N');
        v_('N-4') = -4+v_('N');
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
        for I=v_('1'):v_('N')
            v_('I+1') = 1+I;
            v_('I-5') = -5+I;
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 2.0;
            v_('A') = v_('I-5');
            v_('B') = v_('1');
            v_('A') = v_('A');
            v_('ABSA') = abs(v_('A'));
            v_('ABSA') = fix(v_('ABSA'));
            v_('B') = v_('B');
            v_('ABSB') = abs(v_('B'));
            v_('ABSB') = fix(v_('ABSB'));
            v_('ABSA+B') = v_('ABSA')+v_('ABSB');
            v_('A') = v_('A')+v_('ABSA+B');
            v_('B') = v_('B')+v_('ABSA+B');
            v_('A/B') = fix(v_('A')/v_('B'));
            v_('B/A') = fix(v_('B')/v_('A'));
            v_('SUM') = v_('A/B')+v_('B/A');
            v_('A') = v_('A')*v_('A/B');
            v_('B') = v_('B')*v_('B/A');
            v_('MAXA,B') = v_('A')+v_('B');
            v_('MAXA,B') = fix(v_('MAXA,B')/v_('SUM'));
            v_('MAXA,B') = v_('MAXA,B')-v_('ABSA+B');
            v_('MAXI-5,1') = v_('MAXA,B');
            v_('A') = v_('I+1');
            v_('B') = v_('N');
            v_('A') = -1*v_('A');
            v_('B') = -1*v_('B');
            v_('A') = v_('A');
            v_('ABSA') = abs(v_('A'));
            v_('ABSA') = fix(v_('ABSA'));
            v_('B') = v_('B');
            v_('ABSB') = abs(v_('B'));
            v_('ABSB') = fix(v_('ABSB'));
            v_('ABSA+B') = v_('ABSA')+v_('ABSB');
            v_('A') = v_('A')+v_('ABSA+B');
            v_('B') = v_('B')+v_('ABSA+B');
            v_('A/B') = fix(v_('A')/v_('B'));
            v_('B/A') = fix(v_('B')/v_('A'));
            v_('SUM') = v_('A/B')+v_('B/A');
            v_('A') = v_('A')*v_('A/B');
            v_('B') = v_('B')*v_('B/A');
            v_('MAXA,B') = v_('A')+v_('B');
            v_('MAXA,B') = fix(v_('MAXA,B')/v_('SUM'));
            v_('MINA,B') = v_('ABSA+B')-v_('MAXA,B');
            v_('MINI+1,N') = v_('MINA,B');
            for J=v_('MAXI-5,1'):v_('MINI+1,N')
                [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = 1.0;
            end
        end
        for K=v_('1'):v_('N/2')
            v_('2K') = 2*K;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(K)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2K')))]);
            valA(end+1) = 4.0;
        end
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
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['OBJ',int2str(I)])) = -1.0;
        end
        for K=v_('1'):v_('N/2')
            pbm.gconst(ig_(['C',int2str(K)])) = 3.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 3.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCUBE',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eXEXP',iet_);
        elftv{it}{1} = 'VM';
        elftv{it}{2} = 'VP';
        elftv{it}{3} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['S',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCUBE';
            ielftype(ie) = iet_('eCUBE');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for K=v_('1'):v_('N/2')
            v_('2K') = 2*K;
            v_('2K-1') = -1+v_('2K');
            v_('2K+1') = 1+v_('2K');
            ename = ['P',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXEXP';
            ielftype(ie) = iet_('eXEXP');
            vname = ['X',int2str(round(v_('2K-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VM',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('2K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('2K+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('VP',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL7d3',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            v_('I+1') = 1+I;
            v_('I-5') = -5+I;
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gL7d3';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 5.0;
            v_('A') = v_('I-5');
            v_('B') = v_('1');
            v_('A') = v_('A');
            v_('ABSA') = abs(v_('A'));
            v_('ABSA') = fix(v_('ABSA'));
            v_('B') = v_('B');
            v_('ABSB') = abs(v_('B'));
            v_('ABSB') = fix(v_('ABSB'));
            v_('ABSA+B') = v_('ABSA')+v_('ABSB');
            v_('A') = v_('A')+v_('ABSA+B');
            v_('B') = v_('B')+v_('ABSA+B');
            v_('A/B') = fix(v_('A')/v_('B'));
            v_('B/A') = fix(v_('B')/v_('A'));
            v_('SUM') = v_('A/B')+v_('B/A');
            v_('A') = v_('A')*v_('A/B');
            v_('B') = v_('B')*v_('B/A');
            v_('MAXA,B') = v_('A')+v_('B');
            v_('MAXA,B') = fix(v_('MAXA,B')/v_('SUM'));
            v_('MAXA,B') = v_('MAXA,B')-v_('ABSA+B');
            v_('MAXI-5,1') = v_('MAXA,B');
            v_('A') = v_('I+1');
            v_('B') = v_('N');
            v_('A') = -1*v_('A');
            v_('B') = -1*v_('B');
            v_('A') = v_('A');
            v_('ABSA') = abs(v_('A'));
            v_('ABSA') = fix(v_('ABSA'));
            v_('B') = v_('B');
            v_('ABSB') = abs(v_('B'));
            v_('ABSB') = fix(v_('ABSB'));
            v_('ABSA+B') = v_('ABSA')+v_('ABSB');
            v_('A') = v_('A')+v_('ABSA+B');
            v_('B') = v_('B')+v_('ABSA+B');
            v_('A/B') = fix(v_('A')/v_('B'));
            v_('B/A') = fix(v_('B')/v_('A'));
            v_('SUM') = v_('A/B')+v_('B/A');
            v_('A') = v_('A')*v_('A/B');
            v_('B') = v_('B')*v_('B/A');
            v_('MAXA,B') = v_('A')+v_('B');
            v_('MAXA,B') = fix(v_('MAXA,B')/v_('SUM'));
            v_('MINA,B') = v_('ABSA+B')-v_('MAXA,B');
            v_('MINI+1,N') = v_('MINA,B');
            for J=v_('MAXI-5,1'):v_('MINI+1,N')
                ig = ig_(['OBJ',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.0;
            end
        end
        for K=v_('1'):v_('N/2')
            ig = ig_(['C',int2str(K)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(K)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               6.26382E+04
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AY-V-V';
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
            g_(1,1) = 2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eCUBE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3;
        if(nargout>1)
            g_(1,1) = 3.0*EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eXEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,1) = U_(2,1)+1;
        U_(2,2) = U_(2,2)-1;
        U_(2,3) = U_(2,3)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        EXPW = exp(IV_(2));
        UEXPW = IV_(1)*EXPW;
        varargout{1} = UEXPW;
        if(nargout>1)
            g_(1,1) = EXPW;
            g_(2,1) = UEXPW;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = EXPW;
                H_(2,1) = H_(1,2);
                H_(2,2) = UEXPW;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL7d3'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        Z = abs(GVAR_);
        varargout{1} = Z^(7.0/3.0);
        if(nargout>1)
            g_ = 7.0*sign(GVAR_)*(Z^(4.0/3.0))/3.0;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 28.0*Z^(1.0/3.0)/9.0;
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

