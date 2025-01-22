function varargout = LEVYMONE8(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LEVYMONE8
%    *********
%    A global optimization example due to Levy & Montalvo 
%    This problem is one of the parameterised set LEVYMONT5-LEVYMONT10
% 
%    Source:  problem 8 in
% 
%    A. V. Levy and A. Montalvo
%    "The Tunneling Algorithm for the Global Minimization of Functions"
%    SIAM J. Sci. Stat. Comp. 6(1) 1985 15:29 
%    https://doi.org/10.1137/0906002
% 
%    nonlinear equations version
% 
%    SIF input: Nick Gould, August 2021
% 
%    classification = 'C-CNOR2-AY-5-10'
% 
%    N is the number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   5              $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LEVYMONE8';

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
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('A') = 1.0;
        v_('K') = 10.0;
        v_('L') = 1.0;
        v_('C') = 0.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('RN') = v_('N');
        v_('A-C') = v_('A')-v_('C');
        v_('PI/N') = v_('PI')/v_('RN');
        v_('KPI/N') = v_('K')*v_('PI/N');
        v_('ROOTKPI/N') = sqrt(v_('KPI/N'));
        v_('N/PI') = v_('RN')/v_('PI');
        v_('N/KPI') = v_('N/PI')/v_('K');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['Q',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('L');
            pbm.gscale(ig,1) = v_('N/PI');
            [ig,ig_] = s2mpjlib('ii',['N',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['N',int2str(I)];
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
        for I=v_('1'):v_('N')
            pbm.gconst(ig_(['Q',int2str(I)])) = v_('A-C');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -10.0*ones(pb.n,1);
        pb.xupper = 10.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 8.0*ones(pb.n,1);
        pb.x0(ix_('X1'),1) = -8.0;
        pb.x0(ix_('X2'),1) = 8.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eS2',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'L';
        elftp{it}{2} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'ePS2',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Z';
        elftp{it}{1} = 'L';
        elftp{it}{2} = 'C';
        elftp{it}{3} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eS2';
        ielftype(ie) = iet_('eS2');
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('L',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('L');
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('C',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('C');
        for I=v_('2'):v_('N')
            v_('I-1') = I-v_('1');
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePS2';
            ielftype(ie) = iet_('ePS2');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-10.0,10.0,8.0);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('L',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('L');
            [~,posep] = ismember('C',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('C');
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('A');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['N',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('ROOTKPI/N');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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
        pb.pbclass = 'C-CNOR2-AY-5-10';
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 4.0*atan(1.0e0);
        varargout{1} = pbm;

    case 'eS2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PIL = pbm.efpar(1)*pbm.elpar{iel_}(1);
        V = PIL*EV_(1)+pbm.efpar(1)*pbm.elpar{iel_}(2);
        SINV = sin(V);
        COSV = cos(V);
        varargout{1} = SINV;
        if(nargout>1)
            g_(1,1) = PIL*COSV;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -PIL*PIL*SINV;
                varargout{3} = H_;
            end
        end

    case 'ePS2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PIL = pbm.efpar(1)*pbm.elpar{iel_}(1);
        U = pbm.elpar{iel_}(1)*EV_(2)+pbm.elpar{iel_}(2)-pbm.elpar{iel_}(3);
        V = PIL*EV_(1)+pbm.efpar(1)*pbm.elpar{iel_}(2);
        SINV = sin(V);
        COSV = cos(V);
        varargout{1} = U*SINV;
        if(nargout>1)
            g_(1,1) = PIL*U*COSV;
            g_(2,1) = pbm.elpar{iel_}(1)*SINV;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -PIL*PIL*U*SINV;
                H_(1,2) = pbm.elpar{iel_}(1)*PIL*COSV;
                H_(2,1) = H_(1,2);
                H_(2,2) = 0.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
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

