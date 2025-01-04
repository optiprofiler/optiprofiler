function varargout = CRAGGLVY(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CRAGGLVY
%    *********
%    Extended Cragg and Levy problem.
%    This problem is a sum of m  sets of 5 groups,
%    There are 2m+2 variables. The Hessian matrix is 7-diagonal.
% 
%    Source:  problem 32 in
%    Ph. L. Toint,
%    "Test problems for partially separable optimization and results
%    for the routine PSPMIN",
%    Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
% 
%    See  also Buckley#18
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-COUR2-AY-V-0'
% 
%    M is the number of group sets
% 
%       Alternative values for the SIF file parameters:
% IE M                   1              $-PARAMETER n = 4     original value
% IE M                   4              $-PARAMETER n = 10
% IE M                   24             $-PARAMETER n = 50
% IE M                   49             $-PARAMETER n = 100
% IE M                   249            $-PARAMETER n = 500
% IE M                   499            $-PARAMETER n = 1000
% IE M                   2499           $-PARAMETER n = 5000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CRAGGLVY';

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
            v_('M') = 4;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
        v_('2M') = 2*v_('M');
        v_('N') = 2+v_('2M');
        v_('1') = 1;
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
        for I=v_('1'):v_('M')
            v_('2I') = 2*I;
            v_('2I-1') = -1+v_('2I');
            v_('2I+1') = 1+v_('2I');
            v_('2I+2') = 2+v_('2I');
            [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I')))]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.01;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I+1')))]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I+2')))]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I-1')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I+2')))]);
            valA(end+1) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 2.0*ones(pb.n,1);
        pb.x0(ix_(['X',int2str(round(v_('1')))]),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXPN',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eTANG',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('M')
            v_('2I') = 2*I;
            v_('2I-1') = -1+v_('2I');
            v_('2I+1') = 1+v_('2I');
            v_('2I+2') = 2+v_('2I');
            ename = ['AE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXPN';
            ielftype(ie) = iet_('eEXPN');
            vname = ['X',int2str(round(v_('2I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['CE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eTANG';
            ielftype(ie) = iet_('eTANG');
            vname = ['X',int2str(round(v_('2I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('2I+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],2.0);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        [it,igt_] = s2mpjlib('ii','gL4',igt_);
        [it,igt_] = s2mpjlib('ii','gL6',igt_);
        [it,igt_] = s2mpjlib('ii','gL8',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['A',int2str(I)]);
            pbm.grftype{ig} = 'gL4';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['AE',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['B',int2str(I)]);
            pbm.grftype{ig} = 'gL6';
            ig = ig_(['C',int2str(I)]);
            pbm.grftype{ig} = 'gL4';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CE',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['D',int2str(I)]);
            pbm.grftype{ig} = 'gL8';
            ig = ig_(['F',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(2)            0.0
% LO SOLTN(4)            1.886566
% LO SOLTN(24)           1.5372D+01
% LO SOLTN(29)           3.2270D+01
% LO SOLTN(249)          1.6745D+02
% LO SOLTN(499)          3.3642D+02
% LO SOLTN(2499)         1.6882D+03
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COUR2-AY-V-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

    case 'eEXPN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVAL = exp(EV_(1));
        varargout{1} = FVAL;
        if(nargout>1)
            g_(1,1) = FVAL;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = FVAL;
                varargout{3} = H_;
            end
        end

    case 'eTANG'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        TANU = tan(IV_(1));
        SECU = 1.0/cos(IV_(1));
        SECUSQ = SECU*SECU;
        varargout{1} = TANU;
        if(nargout>1)
            g_(1,1) = SECUSQ;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0*SECUSQ*TANU;
                varargout{3} = U_.'*H_*U_;
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

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
                varargout{3} = H_;
            end
        end

    case 'gL6'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^6;
        if(nargout>1)
            g_ = 6.0*GVAR_^5;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 30.0*GVAR_^4;
                varargout{3} = H_;
            end
        end

    case 'gL8'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^8;
        if(nargout>1)
            g_ = 8.0*GVAR_^7;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 56.0*GVAR_^6;
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

