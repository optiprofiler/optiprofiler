function varargout = DTOC6(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC6
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, 1 control variable and 1 state variable.
% 
%    The problem is convex.
% 
%    Sources: problem 6 in
%    T.F. Coleman and A. Liao,
%    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
%    Control Problems",
%    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    D.M. Murray and S.J. Yakowitz,
%    "The application of optimal contraol methodology to nonlinear programming
%    problems",
%    Mathematical Programming 21, pp. 331-347, 1981.
% 
%    SIF input: Ph. Toint, August 1993
% 
%    classification = 'OOR2-AN-V-V'
% 
%    Problem variants: they are identified by the value of the parameter N.
% 
%    The problem has 2N-1  variables (of which 1 is fixed),
%    and N-1 constraints
% 
%       Alternative values for the SIF file parameters:
% IE N                   11             $-PARAMETER n =   21, m =  10
% IE N                   21             $-PARAMETER n =   41, m =  20
% IE N                   31             $-PARAMETER n =   61, m =  30
% IE N                   41             $-PARAMETER n =   81, m =  40
% IE N                   51             $-PARAMETER n =  101, m =  50
% IE N                   61             $-PARAMETER n =  121, m =  60
% IE N                   71             $-PARAMETER n =  141, m =  70
% IE N                   81             $-PARAMETER n =  161, m =  80
% IE N                   91             $-PARAMETER n =  181, m =  90
% IE N                   101            $-PARAMETER n =  201, m = 100
% IE N                   501            $-PARAMETER n = 1001, m = 500
% IE N                   1001           $-PARAMETER n = 2001, m =1000
% IE N                   5001           $-PARAMETER n =10001, m =5000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DTOC6';

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
            v_('N') = 11;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for T=v_('1'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(T)],ix_);
            pb.xnames{iv} = ['X',int2str(T)];
        end
        for T=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(T)],ix_);
            pb.xnames{iv} = ['Y',int2str(T)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for T=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['OY',int2str(T)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Y',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = 2.0;
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(T)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = 2.0;
        end
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T)];
            iv = ix_(['Y',int2str(round(v_('T+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['Y',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_(['Y',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['Y',int2str(round(v_('1')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXP',iet_);
        elftv{it}{1} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for T=v_('1'):v_('N-1')
            ename = ['E',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
            vname = ['X',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for T=v_('1'):v_('N-1')
            ig = ig_(['OX',int2str(T)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OY',int2str(T)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['TT',int2str(T)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION(  11)      19.80411526774
% LO SOLUTION(  21)      62.49481326823
% LO SOLUTION(  31)      119.0328455446
% LO SOLUTION(  41)      185.8961987565
% LO SOLUTION(  51)      261.1131573312
% LO SOLUTION(  61)      343.3995402405
% LO SOLUTION(  71)      431.8430677467
% LO SOLUTION(  81)      525.7575776948
% LO SOLUTION(  91)      624.6051839906
% LO SOLUTION( 101)      727.9505731659
% LO SOLUTION( 501)      6846.330143698
% LO SOLUTION(1001)      17176.03828316
% LO SOLUTION(5001)      
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-V-V';
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

    case 'eEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EZ = exp(EV_(1));
        varargout{1} = EZ;
        if(nargout>1)
            g_(1,1) = EZ;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = EZ;
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

