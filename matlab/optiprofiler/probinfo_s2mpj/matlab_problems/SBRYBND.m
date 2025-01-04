function varargout = SBRYBND(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SBRYBND
%    *********
%    Broyden banded system of nonlinear equations, considered in the
%    least square sense.
%    NB: scaled version of BRYBND
% 
%    Source: problem 31 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#73 (p. 41) and Toint#18
% 
%    SIF input: Ph. Toint and Nick Gould, Nov 1997.
% 
%    classification = 'C-CSUR2-AN-V-0'
% 
%    N is the number of equations and variables (variable).
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SBRYBND';

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
        v_('ONE') = 1.0;
        v_('KAPPA1') = 2.0;
        v_('KAPPA2') = 5.0;
        v_('KAPPA3') = 1.0;
        v_('LB') = 5;
        v_('UB') = 1;
        v_('RN') = v_('N');
        v_('RN-1') = -1+v_('RN');
        v_('SCAL') = 12.0;
        v_('1') = 1;
        v_('MLB') = -1*v_('LB');
        v_('MUB') = -1*v_('UB');
        v_('LB+1') = 1+v_('LB');
        v_('N-UB') = v_('N')+v_('MUB');
        v_('N-UB-1') = -1+v_('N-UB');
        v_('-KAPPA3') = -1.0*v_('KAPPA3');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('RAT') = v_('RI-1')/v_('RN-1');
            v_('ARG') = v_('RAT')*v_('SCAL');
            v_(['SCALE',int2str(I)]) = exp(v_('ARG'));
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('LB')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('1'):v_('I-1')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('KAP');
            for J=v_('I+1'):v_('I+UB')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
        end
        for I=v_('LB+1'):v_('N-UB-1')
            v_('I-LB') = I+v_('MLB');
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('I-LB'):v_('I-1')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('KAP');
            for J=v_('I+1'):v_('I+UB')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
        end
        for I=v_('N-UB'):v_('N')
            v_('I-LB') = I+v_('MLB');
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('I-LB'):v_('I-1')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('KAP');
            for J=v_('I+1'):v_('N')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('KAP');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('DIV') = v_('ONE')/v_(['SCALE',int2str(I)]);
            pb.x0(ix_(['X',int2str(I)]),1) = v_('DIV');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eCB',iet_);
        elftv{it}{1} = 'V';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            ename = ['E',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['SCALE',int2str(I)]);
            ename = ['Q',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCB';
            ielftype(ie) = iet_('eCB');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['SCALE',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('LB')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('1'):v_('I-1')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('I+UB')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
        end
        for I=v_('LB+1'):v_('N-UB-1')
            v_('I-LB') = I+v_('MLB');
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('I-LB'):v_('I-1')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('I+UB')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
        end
        for I=v_('N-UB'):v_('N')
            v_('I-LB') = I+v_('MLB');
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('I-LB'):v_('I-1')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('N')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSUR2-AN-V-0';
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PP = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        varargout{1} = PP*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = PP*(EV_(1)+EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*PP;
                varargout{3} = H_;
            end
        end

    case 'eCB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        PP = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        varargout{1} = PP*EV_(1)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 3.0*PP*EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 6.0*PP*EV_(1);
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

