function varargout = DIXCHLNG(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DIXCHLNG
%    *********
% 
%    A constrained problem set as a challenge for SQP methods
%    by L.C.W. Dixon at the APMOD91 Conference.
% 
%    Source:
%    L.C.W. Dixon, personnal communication, Jan 1991.
% 
%    SIF input: Ph. Toint, Feb 1991.
% 
%    classification = 'C-CSOR2-AN-10-5'
% 
%    Other parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DIXCHLNG';

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
        v_('1') = 1;
        v_('2') = 2;
        v_('7') = 7;
        v_('9') = 9;
        v_('10') = 10;
        v_('90.0') = 90.0;
        v_('10.1') = 10.1;
        v_('19.8') = 19.8;
        v_('1/90.0') = 1.0/v_('90.0');
        v_('1/10.1') = 1.0/v_('10.1');
        v_('1/19.8') = 1.0/v_('19.8');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('10')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('7')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
            pbm.gscale(ig,1) = 0.01;
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+3')))]);
            valA(end+1) = 1.0;
            pbm.gscale(ig,1) = v_('1/90.0');
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+2')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['E',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
            pbm.gscale(ig,1) = v_('1/10.1');
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+3')))]);
            valA(end+1) = 1.0;
            pbm.gscale(ig,1) = v_('1/10.1');
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('1/19.8');
        end
        for I=v_('2'):v_('2'):v_('10')
            [ig,ig_] = s2mpjlib('ii',['P',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['P',int2str(I)];
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
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0*ones(ngrp,1);
        for I=v_('1'):v_('7')
            pbm.gconst(ig_(['A',int2str(I)])) = 0.0;
            pbm.gconst(ig_(['C',int2str(I)])) = 0.0;
            pbm.gconst(ig_(['G',int2str(I)])) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('X0A') = 2.0;
        v_('X0M') = -1.0;
        for I=v_('1'):v_('2'):v_('9')
            v_('X0') = v_('X0A')*v_('X0M');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('X0');
            v_('1/X0') = 1.0/v_('X0');
            v_('I+1') = 1+I;
            pb.x0(ix_(['X',int2str(round(v_('I+1')))]),1) = v_('1/X0');
            v_('X0A') = 1.0+v_('X0A');
            v_('X0M') = -1.0*v_('X0M');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eS2PR',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'ePR2',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'ePR4',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        [it,iet_] = s2mpjlib( 'ii', 'ePR6',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        [it,iet_] = s2mpjlib( 'ii', 'ePR8',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        [it,iet_] = s2mpjlib( 'ii', 'ePR10',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        elftv{it}{9} = 'V9';
        elftv{it}{10} = 'V10';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('9')
            ename = ['XSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('7')
            v_('I+1') = 1+I;
            v_('I+3') = 3+I;
            ename = ['PR',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eS2PR';
            ielftype(ie) = iet_('eS2PR');
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'PRD2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePR2';
        ielftype(ie) = iet_('ePR2');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PRD4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePR4';
        ielftype(ie) = iet_('ePR4');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PRD6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePR6';
        ielftype(ie) = iet_('ePR6');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PRD8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePR8';
        ielftype(ie) = iet_('ePR8');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PRD10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePR10';
        ielftype(ie) = iet_('ePR10');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('7')
            v_('I+2') = 2+I;
            ig = ig_(['A',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['B',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['C',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('I+2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['D',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['E',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['F',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PR',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for I=v_('2'):v_('2'):v_('10')
            ig = ig_(['P',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PRD',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CSOR2-AN-10-5';
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

    case 'eS2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-1.0)*(EV_(2)-1.0);
        if(nargout>1)
            g_(1,1) = EV_(2)-1.0;
            g_(2,1) = EV_(1)-1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'ePR2'

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

    case 'ePR4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*EV_(4);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*EV_(4);
            g_(2,1) = EV_(1)*EV_(3)*EV_(4);
            g_(3,1) = EV_(1)*EV_(2)*EV_(4);
            g_(4,1) = EV_(1)*EV_(2)*EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = EV_(3)*EV_(4);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4);
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3);
                H_(4,1) = H_(1,4);
                H_(2,3) = EV_(1)*EV_(4);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3);
                H_(4,2) = H_(2,4);
                H_(3,4) = EV_(1)*EV_(2);
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    case 'ePR6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6);
            g_(2,1) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6);
            g_(3,1) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6);
            g_(4,1) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6);
            g_(5,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6);
            g_(6,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(6,6);
                H_(1,2) = EV_(3)*EV_(4)*EV_(5)*EV_(6);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4)*EV_(5)*EV_(6);
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3)*EV_(5)*EV_(6);
                H_(4,1) = H_(1,4);
                H_(1,5) = EV_(2)*EV_(3)*EV_(4)*EV_(6);
                H_(5,1) = H_(1,5);
                H_(1,6) = EV_(2)*EV_(3)*EV_(4)*EV_(5);
                H_(6,1) = H_(1,6);
                H_(2,3) = EV_(1)*EV_(4)*EV_(5)*EV_(6);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3)*EV_(5)*EV_(6);
                H_(4,2) = H_(2,4);
                H_(2,5) = EV_(1)*EV_(3)*EV_(4)*EV_(6);
                H_(5,2) = H_(2,5);
                H_(2,6) = EV_(1)*EV_(3)*EV_(4)*EV_(5);
                H_(6,2) = H_(2,6);
                H_(3,4) = EV_(1)*EV_(2)*EV_(5)*EV_(6);
                H_(4,3) = H_(3,4);
                H_(3,5) = EV_(1)*EV_(2)*EV_(4)*EV_(6);
                H_(5,3) = H_(3,5);
                H_(3,6) = EV_(1)*EV_(2)*EV_(4)*EV_(5);
                H_(6,3) = H_(3,6);
                H_(4,5) = EV_(1)*EV_(2)*EV_(3)*EV_(6);
                H_(5,4) = H_(4,5);
                H_(4,6) = EV_(1)*EV_(2)*EV_(3)*EV_(5);
                H_(6,4) = H_(4,6);
                H_(5,6) = EV_(1)*EV_(2)*EV_(3)*EV_(4);
                H_(6,5) = H_(5,6);
                varargout{3} = H_;
            end
        end

    case 'ePR8'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
            g_(2,1) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
            g_(3,1) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
            g_(4,1) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
            g_(5,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8);
            g_(6,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8);
            g_(7,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8);
            g_(8,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(8,8);
                H_(1,2) = EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(4,1) = H_(1,4);
                H_(1,5) = EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8);
                H_(5,1) = H_(1,5);
                H_(1,6) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8);
                H_(6,1) = H_(1,6);
                H_(1,7) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8);
                H_(7,1) = H_(1,7);
                H_(1,8) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7);
                H_(8,1) = H_(1,8);
                H_(2,3) = EV_(1)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(4,2) = H_(2,4);
                H_(2,5) = EV_(1)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8);
                H_(5,2) = H_(2,5);
                H_(2,6) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8);
                H_(6,2) = H_(2,6);
                H_(2,7) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8);
                H_(7,2) = H_(2,7);
                H_(2,8) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7);
                H_(8,2) = H_(2,8);
                H_(3,4) = EV_(1)*EV_(2)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(4,3) = H_(3,4);
                H_(3,5) = EV_(1)*EV_(2)*EV_(4)*EV_(6)*EV_(7)*EV_(8);
                H_(5,3) = H_(3,5);
                H_(3,6) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(7)*EV_(8);
                H_(6,3) = H_(3,6);
                H_(3,7) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(8);
                H_(7,3) = H_(3,7);
                H_(3,8) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7);
                H_(8,3) = H_(3,8);
                H_(4,5) = EV_(1)*EV_(2)*EV_(3)*EV_(6)*EV_(7)*EV_(8);
                H_(5,4) = H_(4,5);
                H_(4,6) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(7)*EV_(8);
                H_(6,4) = H_(4,6);
                H_(4,7) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(8);
                H_(7,4) = H_(4,7);
                H_(4,8) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7);
                H_(8,4) = H_(4,8);
                H_(5,6) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(7)*EV_(8);
                H_(6,5) = H_(5,6);
                H_(5,7) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(8);
                H_(7,5) = H_(5,7);
                H_(5,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7);
                H_(8,5) = H_(5,8);
                H_(6,7) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(8);
                H_(7,6) = H_(6,7);
                H_(6,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7);
                H_(8,6) = H_(6,8);
                H_(7,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6);
                H_(8,7) = H_(7,8);
                varargout{3} = H_;
            end
        end

    case 'ePR10'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*...
             EV_(9)*EV_(10);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(2,1) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(3,1) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(4,1) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(5,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(6,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
            g_(7,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
            g_(8,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
            g_(9,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
            g_(10,1) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(10,10);
                H_(1,2) = EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(4,1) = H_(1,4);
                H_(1,5) = EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(5,1) = H_(1,5);
                H_(1,6) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(6,1) = H_(1,6);
                H_(1,7) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
                H_(7,1) = H_(1,7);
                H_(1,8) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
                H_(8,1) = H_(1,8);
                H_(1,9) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
                H_(9,1) = H_(1,9);
                H_(1,10) = EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
                H_(10,1) = H_(1,10);
                H_(2,3) = EV_(1)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(4,2) = H_(2,4);
                H_(2,5) = EV_(1)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(5,2) = H_(2,5);
                H_(2,6) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(6,2) = H_(2,6);
                H_(2,7) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
                H_(7,2) = H_(2,7);
                H_(2,8) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
                H_(8,2) = H_(2,8);
                H_(2,9) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
                H_(9,2) = H_(2,9);
                H_(2,10) = EV_(1)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
                H_(10,2) = H_(2,10);
                H_(3,4) = EV_(1)*EV_(2)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(4,3) = H_(3,4);
                H_(3,5) = EV_(1)*EV_(2)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(5,3) = H_(3,5);
                H_(3,6) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(6,3) = H_(3,6);
                H_(3,7) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
                H_(7,3) = H_(3,7);
                H_(3,8) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
                H_(8,3) = H_(3,8);
                H_(3,9) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
                H_(9,3) = H_(3,9);
                H_(3,10) = EV_(1)*EV_(2)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
                H_(10,3) = H_(3,10);
                H_(4,5) = EV_(1)*EV_(2)*EV_(3)*EV_(6)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(5,4) = H_(4,5);
                H_(4,6) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(6,4) = H_(4,6);
                H_(4,7) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
                H_(7,4) = H_(4,7);
                H_(4,8) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
                H_(8,4) = H_(4,8);
                H_(4,9) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
                H_(9,4) = H_(4,9);
                H_(4,10) = EV_(1)*EV_(2)*EV_(3)*EV_(5)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
                H_(10,4) = H_(4,10);
                H_(5,6) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(7)*EV_(8)*EV_(9)*EV_(10);
                H_(6,5) = H_(5,6);
                H_(5,7) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(8)*EV_(9)*EV_(10);
                H_(7,5) = H_(5,7);
                H_(5,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(9)*EV_(10);
                H_(8,5) = H_(5,8);
                H_(5,9) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(10);
                H_(9,5) = H_(5,9);
                H_(5,10) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(6)*EV_(7)*EV_(8)*EV_(9);
                H_(10,5) = H_(5,10);
                H_(6,7) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(8)*EV_(9)*EV_(10);
                H_(7,6) = H_(6,7);
                H_(6,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(9)*EV_(10);
                H_(8,6) = H_(6,8);
                H_(6,9) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(10);
                H_(9,6) = H_(6,9);
                H_(6,10) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(7)*EV_(8)*EV_(9);
                H_(10,6) = H_(6,10);
                H_(7,8) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(9)*EV_(10);
                H_(8,7) = H_(7,8);
                H_(7,9) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(10);
                H_(9,7) = H_(7,9);
                H_(7,10) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(8)*EV_(9);
                H_(10,7) = H_(7,10);
                H_(8,9) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(10);
                H_(9,8) = H_(8,9);
                H_(8,10) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(9);
                H_(10,8) = H_(8,10);
                H_(9,10) = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5)*EV_(6)*EV_(7)*EV_(8);
                H_(10,9) = H_(9,10);
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

