function varargout = HS93(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS93
%    *********
% 
%    A transformer design problem.
% 
%    Source: problem 93 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'C-COOR2-MY-6-2'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS93';

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
        v_('N') = 6;
        v_('1') = 1;
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
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C2';
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
        pbm.gconst(ig_('C1')) = 2.07e+0;
        pbm.gconst(ig_('C2')) = 1.0e+0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 5.54;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 5.54;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 4.4;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 4.4;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 12.02;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 12.02;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 11.82;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 11.82;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 0.702;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 0.702;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 0.852;
        else
            pb.y0(find(pbm.congrps==ig('X6')),1) = 0.852;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOE1',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        [it,iet_] = s2mpjlib( 'ii', 'eOE2',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        [it,iet_] = s2mpjlib( 'ii', 'eOE3',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        [it,iet_] = s2mpjlib( 'ii', 'eOE4',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X6';
        [it,iet_] = s2mpjlib( 'ii', 'eC1E1',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftv{it}{6} = 'X6';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'OE1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE1';
        ielftype(ie) = iet_('eOE1');
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
        ename = 'OE2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE2';
        ielftype(ie) = iet_('eOE2');
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
        ename = 'OE3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE3';
        ielftype(ie) = iet_('eOE3');
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
        ename = 'OE4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eOE4';
        ielftype(ie) = iet_('eOE4');
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
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C1E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eC1E1';
        ielftype(ie) = iet_('eC1E1');
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
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.04e-2;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OE2');
        pbm.grelw{ig}(posel) = 1.87e-2;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.07e-2;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OE4');
        pbm.grelw{ig}(posel) = 4.37e-2;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0e-3;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 6.2e-4;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OE4');
        pbm.grelw{ig}(posel) = 5.8e-4;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               135.075961
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-MY-6-2';
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

    case 'eOE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(1,1) = U_(1,1)+1;
        U_(2,4) = U_(2,4)+1;
        U_(3,1) = U_(3,1)+1;
        U_(3,2) = U_(3,2)+1;
        U_(3,3) = U_(3,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*IV_(3);
        if(nargout>1)
            g_(1,1) = IV_(2)*IV_(3);
            g_(2,1) = IV_(1)*IV_(3);
            g_(3,1) = IV_(1)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = IV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = IV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eOE2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(1,2) = U_(1,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(3,1) = U_(3,1)+1;
        U_(3,2) = U_(3,2)+1.570000e+00;
        U_(3,4) = U_(3,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*IV_(3);
        if(nargout>1)
            g_(1,1) = IV_(2)*IV_(3);
            g_(2,1) = IV_(1)*IV_(3);
            g_(3,1) = IV_(1)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = IV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = IV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eOE3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(4,5);
        U_(1,1) = U_(1,1)+1;
        U_(2,4) = U_(2,4)+1;
        U_(3,5) = U_(3,5)+1;
        U_(4,1) = U_(4,1)+1;
        U_(4,2) = U_(4,2)+1;
        U_(4,3) = U_(4,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        IV_(4) = U_(4,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*(IV_(3)^2)*IV_(4);
        if(nargout>1)
            g_(1,1) = IV_(2)*(IV_(3)^2)*IV_(4);
            g_(2,1) = IV_(1)*(IV_(3)^2)*IV_(4);
            g_(3,1) = IV_(1)*IV_(2)*2.0e+0*IV_(3)*IV_(4);
            g_(4,1) = IV_(1)*IV_(2)*(IV_(3)^2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = (IV_(3)^2)*IV_(4);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2)*2.0e+0*IV_(3)*IV_(4);
                H_(3,1) = H_(1,3);
                H_(1,4) = IV_(2)*(IV_(3)^2);
                H_(4,1) = H_(1,4);
                H_(2,3) = IV_(1)*2.0e+0*IV_(3)*IV_(4);
                H_(3,2) = H_(2,3);
                H_(2,4) = IV_(1)*(IV_(3)^2);
                H_(4,2) = H_(2,4);
                H_(3,3) = IV_(1)*IV_(2)*2.0e+0*IV_(4);
                H_(3,4) = IV_(1)*IV_(2)*2.0e+0*IV_(3);
                H_(4,3) = H_(3,4);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eOE4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(4,5);
        U_(1,2) = U_(1,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(3,5) = U_(3,5)+1;
        U_(4,1) = U_(4,1)+1;
        U_(4,2) = U_(4,2)+1.570000e+00;
        U_(4,4) = U_(4,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        IV_(4) = U_(4,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*(IV_(3)^2)*IV_(4);
        if(nargout>1)
            g_(1,1) = IV_(2)*(IV_(3)^2)*IV_(4);
            g_(2,1) = IV_(1)*(IV_(3)^2)*IV_(4);
            g_(3,1) = IV_(1)*IV_(2)*2.0e+0*IV_(3)*IV_(4);
            g_(4,1) = IV_(1)*IV_(2)*(IV_(3)^2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = (IV_(3)^2)*IV_(4);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2)*2.0e+0*IV_(3)*IV_(4);
                H_(3,1) = H_(1,3);
                H_(1,4) = IV_(2)*(IV_(3)^2);
                H_(4,1) = H_(1,4);
                H_(2,3) = IV_(1)*2.0e+0*IV_(3)*IV_(4);
                H_(3,2) = H_(2,3);
                H_(2,4) = IV_(1)*(IV_(3)^2);
                H_(4,2) = H_(2,4);
                H_(3,3) = IV_(1)*IV_(2)*2.0e+0*IV_(4);
                H_(3,4) = IV_(1)*IV_(2)*2.0e+0*IV_(3);
                H_(4,3) = H_(3,4);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eC1E1'

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

