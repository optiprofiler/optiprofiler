function varargout = HS87(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS87
%    *********
% 
%    Optimization of an electrical network (EDF) by P. Huard.
% 
%    Source: problem 87 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    Note: There are two variants described in the papers
% 
%       D.H. Himmelblau "Applied nonlinear programming",
%       McGraw-Hill, New-York, 1972, problem 15,
% 
%    and
% 
%       A.R. Colville, "A comparative study on nonlinear programming",
%       IBM Scientific Center Report 320-2949, New York, 1968, problem 6.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'C-COOI2-MN-6-4'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS87';

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
        v_('A') = 131.078;
        v_('B') = 1.48577;
        v_('C') = 0.90798;
        v_('D') = cos(1.47588);
        v_('E') = sin(1.47588);
        v_('F') = 1.48577;
        v_('1') = 1;
        v_('-B') = -1.0e+0*v_('B');
        v_('-F') = -1.0e+0*v_('F');
        v_('C/A') = v_('C')/v_('A');
        v_('1/A') = 1.0e+0/v_('A');
        v_('-1/A') = -1.0e+0*v_('1/A');
        v_('CD/A') = v_('C/A')*v_('D');
        v_('CE/A') = v_('C/A')*v_('E');
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
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0e+0;
        end
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C4';
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
        pbm.gconst(ig_('C1')) = -3.0e+2;
        pbm.gconst(ig_('C4')) = -2.0e+2;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X3'),1) = 340.0;
        pb.xlower(ix_('X4'),1) = 340.0;
        pb.xlower(ix_('X5'),1) = -1000.0;
        pb.xupper(ix_('X1')) = 400.0;
        pb.xupper(ix_('X2')) = 1000.0;
        pb.xupper(ix_('X3')) = 420.0;
        pb.xupper(ix_('X4')) = 420.0;
        pb.xupper(ix_('X5')) = 10000.0;
        pb.xupper(ix_('X6')) = 0.5236;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 107.8119;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 107.8119;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 196.3186;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 196.3186;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 373.8307;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 373.8307;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 420.0;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 420.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 21.30713;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 21.30713;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 0.153292;
        else
            pb.y0(find(pbm.congrps==ig('X6')),1) = 0.153292;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eF1',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eF2',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCOS',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eSIN',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eSQUARE',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OF1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eF1';
        ielftype(ie) = iet_('eF1');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OF2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eF2';
        ielftype(ie) = iet_('eF2');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C1E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOS';
        ielftype(ie) = iet_('eCOS');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('-F');
        ename = 'C1E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQUARE';
        ielftype(ie) = iet_('eSQUARE');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C2E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOS';
        ielftype(ie) = iet_('eCOS');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('B');
        ename = 'C2E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQUARE';
        ielftype(ie) = iet_('eSQUARE');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C3E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSIN';
        ielftype(ie) = iet_('eSIN');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('B');
        ename = 'C3E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQUARE';
        ielftype(ie) = iet_('eSQUARE');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'C4E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSIN';
        ielftype(ie) = iet_('eSIN');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('-B');
        ename = 'C4E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQUARE';
        ielftype(ie) = iet_('eSQUARE');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OF1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('OF2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('CD/A');
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C2E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C2E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('CD/A');
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C3E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-1/A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C3E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('CE/A');
        ig = ig_('C4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C4E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('1/A');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C4E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('CE/A');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               8927.5977
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOI2-MN-6-4';
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

    case 'eF1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        I1 = EV_(1)<300.0;
        I2 = EV_(1)>=300.0;
        if(I1)
            F = 30.0*EV_(1);
        end
        if(I2)
            F = 31.0*EV_(1);
        end
        if(I1)
            G = 30.0;
        end
        if(I2)
            G = 31.0;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 0.0e+0;
                varargout{3} = H_;
            end
        end

    case 'eF2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        I1 = EV_(1)<100.0;
        I2 = EV_(1)>=100.0&&EV_(1)<200.0;
        I3 = EV_(1)>=200.0;
        if(I1)
            F = 28.0*EV_(1);
        end
        if(I2)
            F = 29.0*EV_(1);
        end
        if(I3)
            F = 30.0*EV_(1);
        end
        if(I1)
            G = 28.0;
        end
        if(I2)
            G = 29.0;
        end
        if(I3)
            G = 30.0;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 0.0e+0;
                varargout{3} = H_;
            end
        end

    case 'eCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SN = sin(EV_(3)+pbm.elpar{iel_}(1));
        CS = cos(EV_(3)+pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*EV_(2)*CS;
        if(nargout>1)
            g_(1,1) = EV_(2)*CS;
            g_(2,1) = EV_(1)*CS;
            g_(3,1) = -EV_(1)*EV_(2)*SN;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = CS;
                H_(2,1) = H_(1,2);
                H_(1,3) = -EV_(2)*SN;
                H_(3,1) = H_(1,3);
                H_(2,3) = -EV_(1)*SN;
                H_(3,2) = H_(2,3);
                H_(3,3) = -EV_(1)*EV_(2)*CS;
                varargout{3} = H_;
            end
        end

    case 'eSIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SN = sin(EV_(3)+pbm.elpar{iel_}(1));
        CS = cos(EV_(3)+pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*EV_(2)*SN;
        if(nargout>1)
            g_(1,1) = EV_(2)*SN;
            g_(2,1) = EV_(1)*SN;
            g_(3,1) = EV_(1)*EV_(2)*CS;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = SN;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*CS;
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1)*CS;
                H_(3,2) = H_(2,3);
                H_(3,3) = -EV_(1)*EV_(2)*SN;
                varargout{3} = H_;
            end
        end

    case 'eSQUARE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0;
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

