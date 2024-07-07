function varargout = HS81(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS81
%    *******
% 
%    Source: problem 80 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: B.Debarsy, Jan 1990.
% 
%    classification = 'OOR2-AN-5-3'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS81';

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
        v_('N') = 5;
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
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
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
        pbm.gconst(ig_('OBJ')) = 0.5;
        pbm.gconst(ig_('C1')) = 10.0;
        pbm.gconst(ig_('C3')) = -1.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = -2.3;
        pb.xupper(ix_('X1')) = 2.3;
        pb.xlower(ix_('X2'),1) = -2.3;
        pb.xupper(ix_('X2')) = 2.3;
        pb.xlower(ix_('X3'),1) = -3.2;
        pb.xupper(ix_('X3')) = 3.2;
        pb.xlower(ix_('X4'),1) = -3.2;
        pb.xupper(ix_('X4')) = 3.2;
        pb.xlower(ix_('X5'),1) = -3.2;
        pb.xupper(ix_('X5')) = 3.2;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = -2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = -2.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 2.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 2.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = -1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = -1.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = -1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXPO',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        [it,iet_] = s2mpjlib( 'ii', 'eX6',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eX3Y3',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eCB',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'OBJ1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEXPO';
        ielftype(ie) = iet_('eEXPO');
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
        ename = 'OBJ2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX6';
        ielftype(ie) = iet_('eX6');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OBJ3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX6';
        ielftype(ie) = iet_('eX6');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OBJ4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eX3Y3';
        ielftype(ie) = iet_('eX3Y3');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OBJ5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCB';
        ielftype(ie) = iet_('eCB');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OBJ6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCB';
        ielftype(ie) = iet_('eCB');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('1'):v_('N')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCB';
        ielftype(ie) = iet_('eCB');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCB';
        ielftype(ie) = iet_('eCB');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OBJ1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OBJ2');
        pbm.grelw{ig}(posel) = -0.5;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OBJ3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OBJ4');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OBJ5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OBJ6');
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        pbm.grelw{ig}(posel) = -5.0;
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.539498
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-5-3';
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

    case 'eEXPO'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PROD = EV_(1)*EV_(2)*EV_(3)*EV_(4)*EV_(5);
        EX = exp(PROD);
        EX1 = EX*(1.0+PROD);
        varargout{1} = EX;
        if(nargout>1)
            g_(1,1) = EX*EV_(2)*EV_(3)*EV_(4)*EV_(5);
            g_(2,1) = EX*EV_(1)*EV_(3)*EV_(4)*EV_(5);
            g_(3,1) = EX*EV_(1)*EV_(2)*EV_(4)*EV_(5);
            g_(4,1) = EX*EV_(1)*EV_(2)*EV_(3)*EV_(5);
            g_(5,1) = EX*EV_(1)*EV_(2)*EV_(3)*EV_(4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = EX*(EV_(2)*EV_(3)*EV_(4)*EV_(5))^2;
                H_(1,2) = EX1*EV_(3)*EV_(4)*EV_(5);
                H_(2,1) = H_(1,2);
                H_(1,3) = EX1*EV_(2)*EV_(4)*EV_(5);
                H_(3,1) = H_(1,3);
                H_(1,4) = EX1*EV_(2)*EV_(3)*EV_(5);
                H_(4,1) = H_(1,4);
                H_(1,5) = EX1*EV_(2)*EV_(3)*EV_(4);
                H_(5,1) = H_(1,5);
                H_(2,2) = EX*(EV_(1)*EV_(3)*EV_(4)*EV_(5))^2;
                H_(2,3) = EX1*EV_(1)*EV_(4)*EV_(5);
                H_(3,2) = H_(2,3);
                H_(2,4) = EX1*EV_(1)*EV_(3)*EV_(5);
                H_(4,2) = H_(2,4);
                H_(2,5) = EX1*EV_(1)*EV_(3)*EV_(4);
                H_(5,2) = H_(2,5);
                H_(3,3) = EX*(EV_(1)*EV_(2)*EV_(4)*EV_(5))^2;
                H_(3,4) = EX1*EV_(1)*EV_(2)*EV_(5);
                H_(4,3) = H_(3,4);
                H_(3,5) = EX1*EV_(1)*EV_(2)*EV_(4);
                H_(5,3) = H_(3,5);
                H_(4,4) = EX*(EV_(1)*EV_(2)*EV_(3)*EV_(5))^2;
                H_(4,5) = EX1*EV_(1)*EV_(2)*EV_(3);
                H_(5,4) = H_(4,5);
                H_(5,5) = EX*(EV_(1)*EV_(2)*EV_(3)*EV_(4))^2;
                varargout{3} = H_;
            end
        end

    case 'eX6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^6;
        if(nargout>1)
            g_(1,1) = 6.0*EV_(1)^5;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 30.0*EV_(1)^4;
                varargout{3} = H_;
            end
        end

    case 'eX3Y3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        X2 = EV_(1)*EV_(1);
        X3 = X2*EV_(1);
        Y2 = EV_(2)*EV_(2);
        Y3 = Y2*EV_(2);
        varargout{1} = X3*Y3;
        if(nargout>1)
            g_(1,1) = 3.0*X2*Y3;
            g_(2,1) = 3.0*X3*Y2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 6.0*EV_(1)*Y3;
                H_(1,2) = 9.0*X2*Y2;
                H_(2,1) = H_(1,2);
                H_(2,2) = 6.0*X3*EV_(2);
                varargout{3} = H_;
            end
        end

    case 'eCB'

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

