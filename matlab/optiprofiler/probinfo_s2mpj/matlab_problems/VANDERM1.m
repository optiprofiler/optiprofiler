function varargout = VANDERM1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A nonlinear equation problem, subject to monotonicity constraints.
%    The Jacobian is a dense Vandermonde matrix.
% 
%    Problems VANDERM1, VANDERM2, VANDERM3 and VANDERM4 differ by the rhs
%    of the equation.  They are increasingly degenerate.
% 
%    The problem is non-convex.
% 
%    Source:
%    A. Neumaier, private communication, 1991.
% 
%    SIF input: Ph. L. Toint, May 1993.
% 
%    classification = 'NOR2-AN-V-V'
% 
%    Size of the system
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER
% IE N                   3              $-PARAMETER
% IE N                   4              $-PARAMETER
% IE N                   5              $-PARAMETER
% IE N                   10             $-PARAMETER     original value
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'VANDERM1';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_(['AL',int2str(I)]) = v_('RI')/v_('RN');
        end
        v_(['A',int2str(round(v_('1')))]) = 0.0;
        for I=v_('1'):v_('N')
            v_(['A',int2str(round(v_('1')))]) =...
                  v_(['AL',int2str(I)])+v_(['A',int2str(round(v_('1')))]);
        end
        for K=v_('2'):v_('N')
            v_('RK') = K;
            v_(['A',int2str(K)]) = 0.0;
            for I=v_('1'):v_('N')
                v_('LOGAL') = log(v_(['AL',int2str(I)]));
                v_('KLOGAL') = v_('LOGAL')*v_('RK');
                v_('ALK') = exp(v_('KLOGAL'));
                v_(['A',int2str(K)]) = v_(['A',int2str(K)])+v_('ALK');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('1')))];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for K=v_('2'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(K)];
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['M',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['M',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        for K=v_('1'):v_('N')
            pbm.gconst(ig_(['E',int2str(K)])) = v_(['A',int2str(K)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TMP') = v_('RI-1')/v_('RN');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('TMP');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePOWER',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'PWR';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for K=v_('2'):v_('N')
            v_('RK') = K;
            for I=v_('1'):v_('N')
                ename = ['E',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePOWER';
                ielftype(ie) = iet_('ePOWER');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('PWR',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RK');
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_(['E',int2str(round(v_('1')))]);
        pbm.grftype{ig} = 'gL2';
        for K=v_('2'):v_('N')
            ig = ig_(['E',int2str(K)]);
            pbm.grftype{ig} = 'gL2';
            for I=v_('1'):v_('N')
                ig = ig_(['E',int2str(K)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION            0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePOWER'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^pbm.elpar{iel_}(1);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(1)^(pbm.elpar{iel_}(1)-1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) =...
                      pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0)*EV_(1)^(pbm.elpar{iel_}(1)-2.0);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

