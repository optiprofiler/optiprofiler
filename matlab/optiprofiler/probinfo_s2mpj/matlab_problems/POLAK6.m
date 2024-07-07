function varargout = POLAK6(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : POLAK6
%    *********
% 
%    A nonlinear minmax problem in four variables. This is a variation
%    on problem ROSENMMX.
% 
%    Source: 
%    E. Polak, D.H. Mayne and J.E. Higgins,
%    "Superlinearly convergent algorithm for min-max problems"
%    JOTA 69, pp. 407-439, 1991.
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'LOR2-AN-5-4'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'POLAK6';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','F1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F1';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -21.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -21.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.0;
        end
        [ig,ig_] = s2mpjlib('ii','F2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F2';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -15.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -11.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -11.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        [ig,ig_] = s2mpjlib('ii','F3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F3';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -15.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -21.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -21.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        [ig,ig_] = s2mpjlib('ii','F4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F4';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 15.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -15.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -15.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -21.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -21.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
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
        pbm.gconst(ig_('F2')) = 80.0;
        pbm.gconst(ig_('F3')) = 100.0;
        pbm.gconst(ig_('F4')) = 50.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'XX';
        [it,iet_] = s2mpjlib( 'ii', 'eEL42',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        [it,iet_] = s2mpjlib( 'ii', 'eEL442',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        elftv{it}{3} = 'ZZ';
        [it,iet_] = s2mpjlib( 'ii', 'eEL4',iet_);
        elftv{it}{1} = 'XX';
        [it,iet_] = s2mpjlib( 'ii', 'eEL44',iet_);
        elftv{it}{1} = 'XX';
        elftv{it}{2} = 'YY';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL42';
        ielftype(ie) = iet_('eEL42');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL442';
        ielftype(ie) = iet_('eEL442');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('ZZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL4';
        ielftype(ie) = iet_('eEL4');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL44';
        ielftype(ie) = iet_('eEL44');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('YY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X3SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'X4SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('F1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X3SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X4SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 5.0;
        ig = ig_('F2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 11.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 11.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X3SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 12.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X4SQ');
        pbm.grelw{ig}(posel) = 11.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -5.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 15.0;
        ig = ig_('F3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 11.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 21.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X3SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 12.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X4SQ');
        pbm.grelw{ig}(posel) = 21.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 15.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 5.0;
        ig = ig_('F4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 11.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 11.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('X3SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 12.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('X4SQ');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -15.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 15.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution at ( 0, 1, 2, -1 )
% LO SOLTN               -44.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-AN-5-4';
        pb.x0          = zeros(pb.n,1);
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

    case 'eEL42'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(2)+1.0;
        A = EV_(1)-B^4;
        varargout{1} = A*A;
        if(nargout>1)
            g_(1,1) = 2.0*A;
            g_(2,1) = -8.0*A*B^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0;
                H_(1,2) = -8.0*B^3;
                H_(2,1) = H_(1,2);
                H_(2,2) = 32.0*B^6-24.0*A*B^2;
                varargout{3} = H_;
            end
        end

    case 'eEL442'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(3)+1.0;
        C = EV_(2)-B^4;
        DCDZ = -4.0*B^3;
        D2CDZZ = -12.0*B^2;
        A = EV_(1)-C^4;
        DADY = -4.0*C^3;
        DADZ = DADY*DCDZ;
        D2ADYY = -12.0*C^2;
        D2ADYZ = D2ADYY*DCDZ;
        D2ADZZ = D2ADYZ*DCDZ+DADY*D2CDZZ;
        varargout{1} = A*A;
        if(nargout>1)
            g_(1,1) = 2.0*A;
            g_(2,1) = 2.0*A*DADY;
            g_(3,1) = 2.0*A*DADZ;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0;
                H_(1,2) = 2.0*DADY;
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0*DADZ;
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*(DADY*DADY+A*D2ADYY);
                H_(2,3) = 2.0*(DADZ*DADY+A*D2ADYZ);
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*(DADZ*DADZ+A*D2ADZZ);
                varargout{3} = H_;
            end
        end

    case 'eEL4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(1)+1.0;
        varargout{1} = B^4;
        if(nargout>1)
            g_(1,1) = 4.0*B^3;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 12.0*B^2;
                varargout{3} = H_;
            end
        end

    case 'eEL44'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(2)+1.0;
        A = EV_(1)-B^4;
        varargout{1} = A^4;
        if(nargout>1)
            g_(1,1) = 4.0*A^3;
            g_(2,1) = -16.0*(A*B)^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 12.0*A^2;
                H_(1,2) = -48.0*A^2*B^3;
                H_(2,1) = H_(1,2);
                H_(2,2) = -48.0*(A*B)^2*(A-4.0*B^4);
                varargout{3} = H_;
            end
        end

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

