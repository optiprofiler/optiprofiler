function varargout = DTOC4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DTOC4
%    *********
% 
%    This is a discrete time optimal control (DTOC) problem.  
%    The system has N time periods, 1 control variable and 2 state variables.
% 
%    The problem is not convex.
% 
%    Sources: problem 4 in
%    T.F. Coleman and A. Liao,
%    "An Efficient Trust Region Method for Unconstrained Discret-Time Optimal
%    Control Problems",
%    Tech. Report, ctc93tr144,  Advanced Computing Research Institute, 
%    Cornell University, 1992.
% 
%    G. Di Pillo, L. Grippo and F. Lampariello,
%    "A class of structures quasi-Newton algorithms for optimal control
%    problems",
%    in H.E. Rauch, ed., IFAC Applications of nonlinear programming to
%    optimization and control, pp. 101-107, IFAC, Pergamon Press, 1983.
% 
%    SIF input: Ph. Toint, August 1993
% 
%    classification = 'QOR2-AN-V-V'
% 
%    Problem variants: they are identified by the value of the parameter N.
% 
%    The problem has 3N-1  variables (of which 2 are fixed),
%    and 2(N-1) constraints
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER  n=   29,m= 18 original value
% IE N                   50             $-PARAMETER  n=  149,m= 98
% IE N                   100            $-PARAMETER  n=  299,m=198
% IE N                   500            $-PARAMETER  n= 1499,m=998
% IE N                   1000           $-PARAMETER  n= 2999,m=1998
% IE N                   1500           $-PARAMETER  n= 4499,m=2998
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DTOC4';

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
% IE N                   5000           $-PARAMETER  n=14999,m=9998
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('5H') = 5.0*v_('H');
        v_('1/5H') = 1.0/v_('5H');
        v_('1+5H') = 1.0+v_('5H');
        v_('-5H') = -1.0*v_('5H');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for T=v_('1'):v_('N-1')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(T)],ix_);
            pb.xnames{iv} = ['X',int2str(T)];
        end
        for T=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(T),',',int2str(round(v_('1')))],ix_);
            pb.xnames{iv} = ['Y',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(T),',',int2str(round(v_('2')))],ix_);
            pb.xnames{iv} = ['Y',int2str(T),',',int2str(round(v_('2')))];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = v_('1/5H');
        for T=v_('1'):v_('N-1')
            v_('T+1') = 1+T;
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1+5H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1+5H');
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-5H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-5H');
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('1')))];
            iv = ix_(['X',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('5H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('5H');
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('2')))];
            iv = ix_(['Y',int2str(round(v_('T+1'))),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['TT',int2str(T),',',int2str(round(v_('2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['TT',int2str(T),',',int2str(round(v_('2')))];
            iv = ix_(['Y',int2str(T),',',int2str(round(v_('1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('5H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('5H');
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
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) = 1.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))]),1) =...
              0.0;
        pb.x0(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))]),1) =...
              1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eAAB',iet_);
        elftv{it}{1} = 'A';
        elftv{it}{2} = 'B';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = ['Y1SQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['Y1SQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('1'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Y2SQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['Y2SQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('1'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['XSQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['XSQ',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for T=v_('2'):v_('N-1')
            ename = ['Y1SQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['Y',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y2SQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['Y',int2str(T),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['XSQ',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(T)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['Y1SQ',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['Y1SQ',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('N'))),',',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['Y2SQ',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        ename = ['Y2SQ',int2str(round(v_('N')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['Y',int2str(round(v_('N'))),',',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for T=v_('1'):v_('N-1')
            ename = ['E',int2str(T)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eAAB';
            ielftype(ie) = iet_('eAAB');
            vname = ['Y',int2str(T),',',int2str(round(v_('2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(T),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Y1SQ',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Y2SQ',int2str(round(v_('1')))]);
        pbm.grelw{ig}(posel) = 0.5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        for T=v_('2'):v_('N-1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Y1SQ',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Y2SQ',int2str(T)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Y1SQ',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Y2SQ',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 0.5;
        for T=v_('1'):v_('N-1')
            ig = ig_(['TT',int2str(T),',',int2str(round(v_('1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(T)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-5H');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
% LO SOLUTION(  10)      3.75078392210
% LO SOLUTION(  50)      3.02963141755
% LO SOLUTION( 100)      2.94726711402
% LO SOLUTION( 500)      2.87827434035
% LO SOLUTION(1000)      2.87483889886
% LO SOLUTION(5000)      2.86386891514
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QOR2-AN-V-V';
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

    case 'eAAB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*EV_(2);
            g_(2,1) = EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = EV_(2)+EV_(2);
                H_(1,2) = EV_(1)+EV_(1);
                H_(2,1) = H_(1,2);
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

