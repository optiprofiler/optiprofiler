function varargout = CLNLBEAM(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CLNLBEAM
%    *********
% 
%    An optimal control version of the CLamped NonLinear BEAM problem.
%    The energy of a beam of length 1 compressed by a force P is to be
%    minimized.  The control variable is the derivative of the deflection angle.
% 
%    The problem is discretized using the trapezoidal rule. It is non-convex.
% 
%    Source:
%    H. Maurer and H.D. Mittelman,
%    "The non-linear beam via optimal control with bound state variables",
%    Optimal Control Applications and Methods 12, pp. 19-31, 1991.
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'C-COOR2-MN-V-V'
% 
%    Discretization: specify the number of interior points + 1
% 
%       Alternative values for the SIF file parameters:
% IE NI                  10             $-PARAMETER n=33, m=20
% IE NI                  50             $-PARAMETER n=153, m=100
% IE NI                  100            $-PARAMETER n=303, m=200
% IE NI                  500            $-PARAMETER n=1503, m=1000
% IE NI                  1000           $-PARAMETER n=3003, m=2000 original value
% IE NI                  2000           $-PARAMETER n=6003, m=4000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CLNLBEAM';

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
        if(nargs<1)
            v_('NI') = 10;  %  SIF file default value
        else
            v_('NI') = varargin{1};
        end
% IE NI                  5000           $-PARAMETER n=15003, m=10000
        if(nargs<2)
            v_('ALPHA') = 350.0;  %  SIF file default value
        else
            v_('ALPHA') = varargin{2};
        end
        v_('RNI') = v_('NI');
        v_('NI-1') = -1+v_('NI');
        v_('H') = 1.0/v_('RNI');
        v_('H/4') = 0.25*v_('H');
        v_('H/2') = 0.5*v_('H');
        v_('AH') = v_('ALPHA')*v_('H');
        v_('AH/2') = 0.5*v_('AH');
        v_('-H/2') = -0.5*v_('H');
        v_('0') = 0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
        end
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','ENERGY',ig_);
        gtype{ig} = '<>';
        for I=v_('0'):v_('NI-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['EX',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EX',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['ET',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ET',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['T',int2str(I)]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('-H/2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = v_('-H/2');
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for I=v_('0'):v_('NI')
            pb.xlower(ix_(['X',int2str(I)]),1) = -0.05;
            pb.xupper(ix_(['X',int2str(I)])) = 0.05;
        end
        for I=v_('0'):v_('NI')
            pb.xlower(ix_(['T',int2str(I)]),1) = -1.0;
            pb.xupper(ix_(['T',int2str(I)])) = 1.0;
        end
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('NI')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('NI')))]),1) = 0.0;
        pb.xlower(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['T',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['T',int2str(round(v_('NI')))]),1) = 0.0;
        pb.xupper(ix_(['T',int2str(round(v_('NI')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('0'):v_('NI')
            v_('RI') = I;
            v_('TT') = v_('RI')*v_('H');
            v_('CTT') = cos(v_('TT'));
            v_('SCTT') = 0.05*v_('CTT');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('SCTT');
            pb.x0(ix_(['T',int2str(I)]),1) = v_('SCTT');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCOS',iet_);
        elftv{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'eSIN',iet_);
        elftv{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('NI')
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOS';
            ielftype(ie) = iet_('eCOS');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['S',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSIN';
            ielftype(ie) = iet_('eSIN');
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['USQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('0'):v_('NI-1')
            v_('I+1') = 1+I;
            ig = ig_(['EX',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-H/2');
            ig = ig_('ENERGY');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['USQ',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['USQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('H/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('AH/2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('AH/2');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           345.0301196587
% LO SOLTN(50)           344.8673691861
% LO SOLTN(100)          344.8801831150
% LO SOLTN(500)          344.8748539754
% LO SOLTN(1000)         344.8788169123
% LO SOLTN(5000)         
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-MN-V-V';
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

    case 'eCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CC = cos(EV_(1));
        SS = sin(EV_(1));
        varargout{1} = CC;
        if(nargout>1)
            g_(1,1) = -SS;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -CC;
                varargout{3} = H_;
            end
        end

    case 'eSIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CC = cos(EV_(1));
        SS = sin(EV_(1));
        varargout{1} = SS;
        if(nargout>1)
            g_(1,1) = CC;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -SS;
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

