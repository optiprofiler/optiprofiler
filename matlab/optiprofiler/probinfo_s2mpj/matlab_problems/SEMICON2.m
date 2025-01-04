function varargout = SEMICON2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SEMICON2
%    *********
% 
%    The semiconductor problem by Rheinboldt, using a finite difference
%    approximation.
% 
%    Source: problem 10 in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-CNOR2-AN-V-V'
% 
%    N  = Number of discretized point inside the interval [a, b]
%    LN = Index of the last negative discretization point
%         (the interest is in the negative part)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SEMICON2';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE LN                  9              $-PARAMETER     original value
        if(nargs<2)
            v_('LN') = 9;  %  SIF file default value
        else
            v_('LN') = varargin{2};
        end
% IE N                   50             $-PARAMETER
% IE LN                  45             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE LN                  90             $-PARAMETER
% IE N                   500            $-PARAMETER
% IE LN                  450            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE LN                  900            $-PARAMETER
% IE N                   5000           $-PARAMETER
% IE LN                  4500           $-PARAMETER
        if(nargs<3)
            v_('LAMBDA') = 0.2;  %  SIF file default value
        else
            v_('LAMBDA') = varargin{3};
        end
        v_('A') = -0.00009;
        v_('B') = 0.00001;
        v_('UA') = 0.0;
        v_('UB') = 700.0;
        v_('CA') = 1.0e12;
        v_('CB') = 1.0e13;
        v_('BETA') = 40.0;
        v_('LN+1') = 1+v_('LN');
        v_('N+1') = 1+v_('N');
        v_('-A') = -1.0*v_('A');
        v_('B-A') = v_('B')+v_('-A');
        v_('RN+1') = v_('N+1');
        v_('TMP') = 1.0/v_('RN+1');
        v_('H') = v_('B-A')*v_('TMP');
        v_('H2') = v_('H')*v_('H');
        v_('LB') = v_('LAMBDA')*v_('BETA');
        v_('H2CA') = v_('H2')*v_('CA');
        v_('H2CB') = v_('H2')*v_('CB');
        v_('LH2CA') = v_('LAMBDA')*v_('H2CA');
        v_('LH2CB') = v_('LAMBDA')*v_('H2CB');
        v_('LUA') = v_('LAMBDA')*v_('UA');
        v_('LUB') = v_('LAMBDA')*v_('UB');
        v_('ULW') = -5.0+v_('LUA');
        v_('UUP') = 5.0+v_('LUB');
        v_('-LB') = -1.0*v_('LB');
        v_('-LUB') = -1.0*v_('LUB');
        v_('-LH2CB') = -1.0*v_('LH2CB');
        v_('0') = 0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('0'):v_('N+1')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(round(v_('I-1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(I)]);
            valA(end+1) = -2.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['U',int2str(round(v_('I+1')))]);
            valA(end+1) = 1.0;
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('LN')
            pbm.gconst(ig_(['G',int2str(I)])) = v_('LH2CA');
        end
        for I=v_('LN+1'):v_('N')
            pbm.gconst(ig_(['G',int2str(I)])) = v_('-LH2CB');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = v_('UUP')*ones(pb.n,1);
        pb.xlower = v_('ULW')*ones(pb.n,1);
        pb.xlower(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.xupper(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.xlower(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        pb.xupper(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        pb.x0(ix_(['U',int2str(round(v_('0')))]),1) = v_('LUA');
        pb.x0(ix_(['U',int2str(round(v_('N+1')))]),1) = v_('LUB');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eWE1',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'LAC';
        elftp{it}{2} = 'LAB';
        elftp{it}{3} = 'LU';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWE1';
            ielftype(ie) = iet_('eWE1');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,v_('ULW'),v_('UUP'),0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('LAC',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LH2CA');
            [~,posep] = ismember('LAB',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-LB');
            [~,posep] = ismember('LU',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LUA');
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWE1';
            ielftype(ie) = iet_('eWE1');
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,v_('ULW'),v_('UUP'),0.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('LAC',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('-LH2CB');
            [~,posep] = ismember('LAB',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LB');
            [~,posep] = ismember('LU',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('LUB');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
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
        pb.pbclass = 'C-CNOR2-AN-V-V';
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

    case 'eWE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVAL =...
              pbm.elpar{iel_}(1)*exp(pbm.elpar{iel_}(2)*(EV_(1)-pbm.elpar{iel_}(3)));
        varargout{1} = FVAL;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(2)*FVAL;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(2)*FVAL;
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

