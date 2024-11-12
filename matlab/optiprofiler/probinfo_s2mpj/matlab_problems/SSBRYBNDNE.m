function varargout = SSBRYBNDNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SSBRYBNDNE
%    *********
%    Broyden banded system of nonlinear equations, considered in the
%    least square sense.
%    NB: scaled version of BRYBND with scaling proposed by Luksan et al.
%    This is a nonlinear equation variant of SSBRYBND
% 
%    Source: problem 48 in
%    L. Luksan, C. Matonoha and J. Vlcek
%    Modified CUTE problems for sparse unconstraoined optimization
%    Technical Report 1081
%    Institute of Computer Science
%    Academy of Science of the Czech Republic
% 
%    that is a scaled variant of problem 31 in
% 
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#73 (p. 41) and Toint#18
% 
%    SIF input: Ph. Toint and Nick Gould, Nov 1997.
%               Nick Gould (nonlinear equation version), Jan 2019
% 
%    classification = 'C-CNOR2-AN-V-V'
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
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SSBRYBNDNE';

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
        v_('SCAL') = 6.0;
        v_('1') = 1;
        v_('MLB') = -1*v_('LB');
        v_('MUB') = -1*v_('UB');
        v_('LB+1') = 1+v_('LB');
        v_('N-UB') = v_('N')+v_('MUB');
        v_('N-UB-1') = -1+v_('N-UB');
        v_('-KAPPA3') = -1.0*v_('KAPPA3');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
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
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('LB')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('1'):v_('I-1')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('KAP');
            end
            for J=v_('I+1'):v_('I+UB')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
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
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('KAP');
            end
            for J=v_('I+1'):v_('I+UB')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
            end
        end
        for I=v_('N-UB'):v_('N')
            v_('I-LB') = I+v_('MLB');
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('I-LB'):v_('I-1')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
            end
            v_('KAP') = v_('KAPPA1')*v_(['SCALE',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('KAP');
            end
            for J=v_('I+1'):v_('N')
                v_('KAP') = v_('-KAPPA3')*v_(['SCALE',int2str(J)]);
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('KAP')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('KAP');
                end
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
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
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
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('LB')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+UB') = I+v_('UB');
            for J=v_('1'):v_('I-1')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('I+UB')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                nlc = union(nlc,ig);
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
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('I+UB')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                nlc = union(nlc,ig);
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
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('KAPPA2');
            for J=v_('I+1'):v_('N')
                ig = ig_(['G',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-KAPPA3');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
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

