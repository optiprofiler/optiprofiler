function varargout = DEVGLA1NE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DEVGLA1NE
%    *********
% 
%    SCIPY global optimization benchmark example DeVilliersGlasser01
% 
%    Fit: y  = x_1 x_2^t sin( t x_3 + x_4 )  +  e
% 
%    Source:  Problem from the SCIPY benchmark set
%      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
%              benchmarks/go_benchmark_functions
% 
%    Nonlinear-equation formulation of DEVGLA1.SIF
% 
%    SIF input: Nick Gould, Jan 2020
% 
%    classification = 'NOR2-MN-4-24'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DEVGLA1NE';

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
        v_('M') = 24;
        v_('N') = 4;
        v_('1') = 1;
        v_('A') = 1.371;
        v_('LNA') = log(v_('A'));
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('RIM1') = -1.0+v_('RI');
            v_('T') = 0.1*v_('RIM1');
            v_(['T',int2str(I)]) = v_('T');
            v_('TLNA') = v_('T')*v_('LNA');
            v_('AT') = exp(v_('TLNA'));
            v_('TP') = 3.112*v_('T');
            v_('TPA') = 1.761+v_('TP');
            v_('STPA') = sin(v_('TPA'));
            v_('P') = v_('AT')*v_('STPA');
            v_('PP') = 60.137*v_('P');
            v_(['Y',int2str(I)]) = v_('PP');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['F',int2str(I)];
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 2.0;
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
            pb.x0(ix_('X4'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eDG1',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDG1';
            ielftype(ie) = iet_('eDG1');
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
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLUTION            0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-4-24';
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

    case 'eDG1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = pbm.elpar{iel_}(1)*EV_(3)+EV_(4);
        SINA = sin(A);
        COSA = cos(A);
        X2T = EV_(2)^pbm.elpar{iel_}(1);
        DX2T = pbm.elpar{iel_}(1)*EV_(2)^(pbm.elpar{iel_}(1)-1.0e0);
        D2X2T =...
              pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0e0)*EV_(2)^(pbm.elpar{iel_}(1)-2.0e0);
        X1X2T = EV_(1)*X2T;
        varargout{1} = X1X2T*SINA;
        if(nargout>1)
            g_(1,1) = X2T*SINA;
            g_(2,1) = EV_(1)*DX2T*SINA;
            g_(3,1) = X1X2T*COSA*pbm.elpar{iel_}(1);
            g_(4,1) = X1X2T*COSA;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = DX2T*SINA;
                H_(2,1) = H_(1,2);
                H_(1,3) = X2T*COSA*pbm.elpar{iel_}(1);
                H_(3,1) = H_(1,3);
                H_(1,4) = X2T*COSA;
                H_(4,1) = H_(1,4);
                H_(2,2) = EV_(1)*D2X2T*SINA;
                H_(2,3) = EV_(1)*DX2T*COSA*pbm.elpar{iel_}(1);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*DX2T*COSA;
                H_(4,2) = H_(2,4);
                H_(3,3) = -X1X2T*SINA*pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
                H_(3,4) = -X1X2T*SINA*pbm.elpar{iel_}(1);
                H_(4,3) = H_(3,4);
                H_(4,4) = -X1X2T*SINA;
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

