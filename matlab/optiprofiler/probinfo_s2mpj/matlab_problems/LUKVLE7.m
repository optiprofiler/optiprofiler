function varargout = LUKVLE7(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKVLE7
%    *********
% 
%    Source: Problem 5.7, A trigonometric tridiagonal function with 
%    simplified five-diagonal constraints, due to L. Luksan and J. Vlcek,
%    "Sparse and partially separable test problems for 
%    unconstrained and equality constrained optimization",
%    Technical Report 767, Inst. Computer Science, Academy of Sciences
%    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
% 
%    SIF input: Nick Gould, April 2001
% 
%    classification = 'OOR2-AY-V-V'
% 
%    some useful parameters, including N, the number of variables.
% 
%       Alternative values for the SIF file parameters:
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   10000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKVLE7';

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
% IE N                   100000         $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('N-3') = -3+v_('N');
        v_('N') = v_('N');
        v_('-N') = -1.0*v_('N');
        v_('N+1') = 1.0+v_('N');
        v_('N.N+1') = v_('N')*v_('N+1');
        v_('-N.N+1/2') = -0.5*v_('N.N+1');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('1')))];
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.0;
        end
        iv = ix_(['X',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('2')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('2')))];
        iv = ix_(['X',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        iv = ix_(['X',int2str(round(v_('3')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('3')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('3')))];
        iv = ix_(['X',int2str(round(v_('N-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        iv = ix_(['X',int2str(round(v_('N-3')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('4')))];
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_(['X',int2str(round(v_('N-2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
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
        pbm.gconst(ig_('OBJ')) = v_('-N.N+1/2');
        pbm.gconst(ig_(['C',int2str(round(v_('2')))])) = 2.0;
        pbm.gconst(ig_(['C',int2str(round(v_('3')))])) = 2.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSIN',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCOS',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCUBEP',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['SI',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSIN';
            ielftype(ie) = iet_('eSIN');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['CO',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOS';
            ielftype(ie) = iet_('eCOS');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['C1',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C1',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['CO',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['SI',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        for I=v_('2'):v_('N-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('I') = I;
            v_('-I') = -1.0*I;
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CO',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-I');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SI',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = I;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SI',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-I');
        end
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['CO',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-N');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['SI',int2str(round(v_('N-1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('N');
        ig = ig_(['C',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_(['C',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               -2.2421E+02
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AY-V-V';
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

    case 'eSIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SINV = sin(EV_(1));
        varargout{1} = SINV;
        if(nargout>1)
            g_(1,1) = cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -SINV;
                varargout{3} = H_;
            end
        end

    case 'eCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        COSV = cos(EV_(1));
        varargout{1} = COSV;
        if(nargout>1)
            g_(1,1) = -sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -COSV;
                varargout{3} = H_;
            end
        end

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

    case 'eCUBEP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3-EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 3.0*EV_(1)^2-EV_(2);
            g_(2,1) = -EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 6.0*EV_(1);
                H_(1,2) = -1.0;
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

