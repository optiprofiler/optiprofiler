function varargout = BRATU2DT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BRATU2DT
%    *********
% 
%    The 2D Bratu problem on the unit square, using finite differences.
%    At the turning point.
% 
%    Source: Problem 3 in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-CNOR2-MN-V-V'
% 
%    P is the number of points in one side of the unit square.
% 
%       Alternative values for the SIF file parameters:
% IE P                   7              $-PARAMETER  n=P**2   original value
% IE P                   10             $-PARAMETER  n=P**2
% IE P                   22             $-PARAMETER  n=P**2
% IE P                   32             $-PARAMETER  n=P**2
% IE P                   72             $-PARAMETER  n=P**2
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BRATU2DT';

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
            v_('P') = 7;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
        if(nargs<2)
            v_('LAMBDA') = 6.80812;  %  SIF file default value
        else
            v_('LAMBDA') = varargin{2};
        end
        v_('1.0') = 1.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('P-1') = -1+v_('P');
        v_('RP-1') = v_('P-1');
        v_('H') = v_('1.0')/v_('RP-1');
        v_('H2') = v_('H')*v_('H');
        v_('C') = v_('H2')*v_('LAMBDA');
        v_('-C') = -1.0*v_('C');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for J=v_('1'):v_('P')
            for I=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('2'):v_('P-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            for J=v_('2'):v_('P-1')
                v_('J+1') = 1+J;
                v_('J-1') = -1+J;
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['U',int2str(I),',',int2str(J)]);
                valA(end+1) = 4.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['U',int2str(round(v_('I+1'))),',',int2str(J)]);
                valA(end+1) = -1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['U',int2str(round(v_('I-1'))),',',int2str(J)]);
                valA(end+1) = -1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['U',int2str(I),',',int2str(round(v_('J+1')))]);
                valA(end+1) = -1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['U',int2str(I),',',int2str(round(v_('J-1')))]);
                valA(end+1) = -1.0;
            end
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
        for J=v_('1'):v_('P')
            pb.xlower(ix_(['U',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['U',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('P-1')
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXP',iet_);
        elftv{it}{1} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eEXP';
                    ielftype(ie) = iet_('eEXP');
                end
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                ig = ig_(['G',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-C');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(4)            1.23159D-02
% LO SOLTN(7)            2.24270D-04
% LO SOLTN(10)           1.85347D-05
% LO SOLTN(22)           1.18376D-07
% LO SOLTN(32)           1.27193D-07
% LO SOLTN(72)           1.30497D-06
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-V-V';
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

    case 'eEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPU = exp(EV_(1));
        varargout{1} = EXPU;
        if(nargout>1)
            g_(1,1) = EXPU;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = EXPU;
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

