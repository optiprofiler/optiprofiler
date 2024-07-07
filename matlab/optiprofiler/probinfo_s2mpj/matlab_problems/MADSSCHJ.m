function varargout = MADSSCHJ(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MADSSCHJ
%    *********
% 
%    A nonlinear minmax problem with variable dimension.
%    The Jacobian of the constraints is dense.
% 
%    Source:
%    K. Madsen and H. Schjaer-Jacobsen,
%    "Linearly Constrained Minmax Optimization",
%    Mathematical Programming 14, pp. 208-223, 1978.
% 
%    SIF input: Ph. Toint, August 1993.
% 
%    classification = 'LQR2-AN-V-V'
% 
%    N is the number of variables - 1, and must be even and at least 4.
%    The number of inequality constraints is 2*N - 2.
% 
%       Alternative values for the SIF file parameters:
% IE N                   4              $-PARAMETER  n=  5, m=  6
% IE N                   10             $-PARAMETER  n= 11, m= 18  original value
% IE N                   20             $-PARAMETER  n= 21, m= 38
% IE N                   30             $-PARAMETER  n= 31, m= 58
% IE N                   40             $-PARAMETER  n= 41, m= 78
% IE N                   50             $-PARAMETER  n= 51, m= 98
% IE N                   60             $-PARAMETER  n= 61, m=118
% IE N                   70             $-PARAMETER  n= 71, m=138
% IE N                   80             $-PARAMETER  n= 81, m=158
% IE N                   90             $-PARAMETER  n= 91, m=178
% IE N                   100            $-PARAMETER  n=101, m=198
% IE N                   200            $-PARAMETER  n=201, m=398
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MADSSCHJ';

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
            v_('N') = 4;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('N-1') = -1+v_('N');
        v_('2N') = v_('N')+v_('N');
        v_('M') = -2+v_('2N');
        v_('M-1') = -1+v_('M');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','Z',ix_);
        pb.xnames{iv} = 'Z';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('2'):v_('N')
            [ig,ig_] = s2mpjlib('ii','C1',ig_);
            gtype{ig}  = '>=';
            cnames{ig} = 'C1';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C2';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('3'):v_('N')
            [ig,ig_] = s2mpjlib('ii','C2',ig_);
            gtype{ig}  = '>=';
            cnames{ig} = 'C2';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C3';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('3'):v_('N')
            [ig,ig_] = s2mpjlib('ii','C3',ig_);
            gtype{ig}  = '>=';
            cnames{ig} = 'C3';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for K=v_('4'):v_('2'):v_('M-1')
            v_('K+1') = 1+K;
            v_('K+2') = 2+K;
            v_('J') = fix(v_('K+2')/v_('2'));
            v_('J-1') = -1+v_('J');
            v_('J+1') = 1+v_('J');
            [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(K)];
            iv = ix_('Z');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('K+1')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(round(v_('K+1')))];
            iv = ix_('Z');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            for I=v_('1'):v_('J-1')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(K)];
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('K+1')))],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(round(v_('K+1')))];
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
            for I=v_('J+1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(K)];
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('K+1')))],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['C',int2str(round(v_('K+1')))];
                iv = ix_(['X',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('M')))],ig_);
        gtype{ig}  = '>=';
        cnames{ig} = ['C',int2str(round(v_('M')))];
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('M')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(round(v_('M')))];
            iv = ix_(['X',int2str(I)]);
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
        for K=v_('1'):v_('M')
            pbm.gconst(ig_(['C',int2str(K)])) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 10.0*ones(pb.n,1);
        pb.x0(ix_('Z'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['XSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],10.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_(['C',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.0;
        for K=v_('4'):v_('2'):v_('M-1')
            v_('K+1') = 1+K;
            v_('K+2') = 2+K;
            v_('J') = fix(v_('K+2')/v_('2'));
            ig = ig_(['C',int2str(K)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('J')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['C',int2str(round(v_('K+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('J')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -2.0;
        end
        ig = ig_(['C',int2str(round(v_('M')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(4)             -2.6121094144
% LO SOLTN(10)            -12.814452425
% LO SOLTN(20)            -49.869888156
% LO SOLTN(30)            -111.93545559
% LO SOLTN(40)            -199.00371592
% LO SOLTN(50)            -311.07308068
% LO SOLTN(60)            -448.14300524
% LO SOLTN(70)            -610.21325256
% LO SOLTN(80)            -797.28370289
% LO SOLTN(90)            -1009.3542892
% LO SOLTN(100)           -1246.4249710
% LO SOLTN(200)           -4992.1339031
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-AN-V-V';
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

