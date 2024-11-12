function varargout = KOWOSB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : KOWOSB
%    *********
% 
%    A problem arising in the analysis of kinetic data for an enzyme
%    reaction, known under the name of Kowalik and Osborne problem
%    in 4 variables.
% 
%    Source:  Problem 15 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'C-CSUR2-MN-4-0'
% 
%    This function  is a nonlinear least squares with 11 groups.  Each
%    group has a linear and a nonlinear element.
% 
%    Number of groups
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'KOWOSB';

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
        v_('M') = 11;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('G1')) = 0.1957;
        pbm.gconst(ig_('G2')) = 0.1947;
        pbm.gconst(ig_('G3')) = 0.1735;
        pbm.gconst(ig_('G4')) = 0.1600;
        pbm.gconst(ig_('G5')) = 0.0844;
        pbm.gconst(ig_('G6')) = 0.0627;
        pbm.gconst(ig_('G7')) = 0.0456;
        pbm.gconst(ig_('G8')) = 0.0342;
        pbm.gconst(ig_('G9')) = 0.0323;
        pbm.gconst(ig_('G10')) = 0.0235;
        pbm.gconst(ig_('G11')) = 0.0246;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 0.25;
        pb.x0(ix_('X2'),1) = 0.39;
        pb.x0(ix_('X3'),1) = 0.415;
        pb.x0(ix_('X4'),1) = 0.39;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eKWO',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftp{it}{1} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eKWO';
            ielftype(ie) = iet_('eKWO');
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
        end
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.0;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.25;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.167;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.125;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.1;
        ename = 'E9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0833;
        ename = 'E10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0714;
        ename = 'E11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        [~,posep] = ismember('U',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0624;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.00102734
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSUR2-MN-4-0';
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

    case 'eKWO'

        EV_  = varargin{1};
        iel_ = varargin{2};
        USQ = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        B1 = USQ+pbm.elpar{iel_}(1)*EV_(2);
        B2 = USQ+pbm.elpar{iel_}(1)*EV_(3)+EV_(4);
        B2SQ = B2*B2;
        B2CB = B2*B2SQ;
        UV1 = pbm.elpar{iel_}(1)*EV_(1);
        UB1 = pbm.elpar{iel_}(1)*B1;
        T1 = B1/B2SQ;
        T2 = 2.0/B2CB;
        varargout{1} = EV_(1)*B1/B2;
        if(nargout>1)
            g_(1,1) = B1/B2;
            g_(2,1) = UV1/B2;
            g_(3,1) = -UV1*T1;
            g_(4,1) = -EV_(1)*T1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = pbm.elpar{iel_}(1)/B2;
                H_(2,1) = H_(1,2);
                H_(1,3) = -UB1/B2SQ;
                H_(3,1) = H_(1,3);
                H_(1,4) = -T1;
                H_(4,1) = H_(1,4);
                H_(2,3) = -UV1*pbm.elpar{iel_}(1)/B2SQ;
                H_(3,2) = H_(2,3);
                H_(2,4) = -UV1/B2SQ;
                H_(4,2) = H_(2,4);
                H_(3,3) = T2*UV1*UB1;
                H_(3,4) = T2*UV1*B1;
                H_(4,3) = H_(3,4);
                H_(4,4) = T2*EV_(1)*B1;
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

