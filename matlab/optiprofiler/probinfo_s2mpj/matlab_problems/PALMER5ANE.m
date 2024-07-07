function varargout = PALMER5ANE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER5ANE
%    *********
% 
%    A nonlinear least squares problem with bounds
%    arising from chemical kinetics.
% 
%    model: H-N=C=Se TZVP + MP2
%    fitting Y to A0 T_0 + A2 T_2 + A4 T_4 + A6 T_6 + A8 T_8 +
%                 A10 T_10 + A12 T_12 + A14 T_14
%                 + B / ( C + X**2 ), B, C nonnegative.
%    where T_i is the i-th (shifted) Chebyshev polynomial
% 
%    Source:
%    M. Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1992.
%    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
% 
%    classification = 'NOR2-RN-8-12'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER5ANE';

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
        v_('M') = 23;
        v_('1') = 1;
        v_('2') = 2;
        v_('12') = 12;
        v_('14') = 14;
        v_('X12') = 0.000000;
        v_('X13') = 1.570796;
        v_('X14') = 1.396263;
        v_('X15') = 1.308997;
        v_('X16') = 1.221730;
        v_('X17') = 1.125835;
        v_('X18') = 1.047198;
        v_('X19') = 0.872665;
        v_('X20') = 0.698132;
        v_('X21') = 0.523599;
        v_('X22') = 0.349066;
        v_('X23') = 0.174533;
        v_('B') = v_('X13');
        v_('A') = -1.0e+0*v_('B');
        v_('DIFF') = 2.0e+0*v_('B');
        v_('Y12') = 83.57418;
        v_('Y13') = 81.007654;
        v_('Y14') = 18.983286;
        v_('Y15') = 8.051067;
        v_('Y16') = 2.044762;
        v_('Y17') = 0.000000;
        v_('Y18') = 1.170451;
        v_('Y19') = 10.479881;
        v_('Y20') = 25.785001;
        v_('Y21') = 44.126844;
        v_('Y22') = 62.822177;
        v_('Y23') = 77.719674;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','A0',ix_);
        pb.xnames{iv} = 'A0';
        [iv,ix_] = s2mpjlib('ii','A2',ix_);
        pb.xnames{iv} = 'A2';
        [iv,ix_] = s2mpjlib('ii','A4',ix_);
        pb.xnames{iv} = 'A4';
        [iv,ix_] = s2mpjlib('ii','A6',ix_);
        pb.xnames{iv} = 'A6';
        [iv,ix_] = s2mpjlib('ii','A8',ix_);
        pb.xnames{iv} = 'A8';
        [iv,ix_] = s2mpjlib('ii','A10',ix_);
        pb.xnames{iv} = 'A10';
        [iv,ix_] = s2mpjlib('ii','B',ix_);
        pb.xnames{iv} = 'B';
        [iv,ix_] = s2mpjlib('ii','C',ix_);
        pb.xnames{iv} = 'C';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('12'):v_('M')
            v_('T0') = 1.0e+0;
            v_('Y') = 2.0e+0*v_(['X',int2str(I)]);
            v_('Y') = v_('Y')-v_('A');
            v_('Y') = v_('Y')-v_('B');
            v_('Y') = v_('Y')/v_('DIFF');
            v_('T1') = v_('Y');
            v_('2Y') = 2.0e+0*v_('Y');
            for J=v_('2'):v_('14')
                v_('J-1') = -1+J;
                v_('J-2') = -2+J;
                v_(['T',int2str(J)]) = v_('2Y')*v_(['T',int2str(round(v_('J-1')))]);
                v_(['T',int2str(J)]) =...
                      v_(['T',int2str(J)])-v_(['T',int2str(round(v_('J-2')))]);
            end
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['O',int2str(I)];
            iv = ix_('A0');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T0')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T0');
            end
            iv = ix_('A2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T2')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T2');
            end
            iv = ix_('A4');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T4')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T4');
            end
            iv = ix_('A6');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T6')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T6');
            end
            iv = ix_('A8');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T8')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T8');
            end
            iv = ix_('A10');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('T10')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('T10');
            end
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
        for I=v_('12'):v_('M')
            pbm.gconst(ig_(['O',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('A0')) = -Inf;
        pb.xupper(ix_('A0'),1) = +Inf;
        pb.xlower(ix_('A2')) = -Inf;
        pb.xupper(ix_('A2'),1) = +Inf;
        pb.xlower(ix_('A4')) = -Inf;
        pb.xupper(ix_('A4'),1) = +Inf;
        pb.xlower(ix_('A6')) = -Inf;
        pb.xupper(ix_('A6'),1) = +Inf;
        pb.xlower(ix_('A8')) = -Inf;
        pb.xupper(ix_('A8'),1) = +Inf;
        pb.xlower(ix_('A10')) = -Inf;
        pb.xupper(ix_('A10'),1) = +Inf;
        pb.xlower(ix_('C'),1) = 0.00001;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eQUOT',iet_);
        elftv{it}{1} = 'B';
        elftv{it}{2} = 'C';
        elftp{it}{1} = 'XSQR';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('12'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eQUOT';
            ielftype(ie) = iet_('eQUOT');
            vname = 'B';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'C';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('C',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('XSQR',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('XSQR');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('12'):v_('M')
            ig = ig_(['O',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
% LO PALMER5A               0.0
%    Solution
% LO SOLTN               4.0606141D-02
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-RN-8-12';
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

    case 'eQUOT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DENOM = 1.0/(EV_(2)+pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*DENOM;
        if(nargout>1)
            g_(1,1) = DENOM;
            g_(2,1) = -EV_(1)*DENOM*DENOM;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -DENOM*DENOM;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1)*DENOM^3;
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

