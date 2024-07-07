function varargout = ORTHREGA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ORTHREGA
%    *********
% 
%    An orthogonal regression problem.
% 
%    The problem is to fit (orthogonally) an ellipse to a set of points
%    in the plane.
% 
%    Source:
%    M. Gulliksson,
%    "Algorithms for nonlinear Least-squares with Applications to
%    Orthogonal Regression",
%    UMINF-178.90, University of Umea, Sweden, 1990.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'QQR2-AN-V-V'
% 
%    Number of levels in the generation of the data points
%    ( number of data points =     4**LEVELS
%      number of variables   = 2 * 4**LEVELS + 5
%      number of constraints =     4**LEVELS         )
% 
%       Alternative values for the SIF file parameters:
% IE LEVELS              3              $-PARAMETER n = 133    original value
% IE LEVELS              4              $-PARAMETER n = 517
% IE LEVELS              5              $-PARAMETER n = 2053
% IE LEVELS              6              $-PARAMETER n = 8197
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORTHREGA';

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
        if(nargs<1)
            v_('LEVELS') = 2;  %  SIF file default value
        else
            v_('LEVELS') = varargin{1};
        end
% IE LEVELS              7              $-PARAMETER n = 32773
% IE LEVELS              8              $-PARAMETER n = 131077
        v_('A') = 9.0;
        v_('B') = 6.0;
        v_('CX') = 0.5;
        v_('CY') = 0.5;
        v_('1') = 1;
        v_('PI') = 3.1415926535;
        v_('-A') = -1.0*v_('A');
        v_('-B') = -1.0*v_('B');
        v_('NPTS') = 1;
        v_(['XD',int2str(round(v_('1')))]) = v_('CX');
        v_(['YD',int2str(round(v_('1')))]) = v_('CY');
        for I=v_('1'):v_('LEVELS')
            v_('NP') = 0+v_('NPTS');
            for J=v_('1'):v_('NP')
                v_(['XZ',int2str(J)]) = v_(['XD',int2str(J)]);
                v_(['YZ',int2str(J)]) = v_(['YD',int2str(J)]);
            end
            v_('NPTS') = 0;
            for J=v_('1'):v_('NP')
                v_('NPTS') = 1+v_('NPTS');
                v_(['XD',int2str(round(v_('NPTS')))]) = v_(['XZ',int2str(J)])+v_('A');
                v_(['YD',int2str(round(v_('NPTS')))]) = v_(['YZ',int2str(J)])+v_('A');
                v_('NPTS') = 1+v_('NPTS');
                v_(['XD',int2str(round(v_('NPTS')))]) = v_(['XZ',int2str(J)])+v_('B');
                v_(['YD',int2str(round(v_('NPTS')))]) = v_(['YZ',int2str(J)])+v_('-B');
                v_('NPTS') = 1+v_('NPTS');
                v_(['XD',int2str(round(v_('NPTS')))]) = v_(['XZ',int2str(J)])+v_('-A');
                v_(['YD',int2str(round(v_('NPTS')))]) = v_(['YZ',int2str(J)])+v_('-A');
                v_('NPTS') = 1+v_('NPTS');
                v_(['XD',int2str(round(v_('NPTS')))]) = v_(['XZ',int2str(J)])+v_('-B');
                v_(['YD',int2str(round(v_('NPTS')))]) = v_(['YZ',int2str(J)])+v_('B');
            end
            v_('A') = v_('A')/v_('PI');
            v_('B') = v_('B')/v_('PI');
            v_('-A') = v_('-A')/v_('PI');
            v_('-B') = v_('-B')/v_('PI');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','H11',ix_);
        pb.xnames{iv} = 'H11';
        [iv,ix_] = s2mpjlib('ii','H12',ix_);
        pb.xnames{iv} = 'H12';
        [iv,ix_] = s2mpjlib('ii','H22',ix_);
        pb.xnames{iv} = 'H22';
        [iv,ix_] = s2mpjlib('ii','G1',ix_);
        pb.xnames{iv} = 'G1';
        [iv,ix_] = s2mpjlib('ii','G2',ix_);
        pb.xnames{iv} = 'G2';
        for I=v_('1'):v_('NPTS')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NPTS')
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['OY',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(I)];
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
        for I=v_('1'):v_('NPTS')
            pbm.gconst(ig_(['OX',int2str(I)])) = v_(['XD',int2str(I)]);
            pbm.gconst(ig_(['OY',int2str(I)])) = v_(['YD',int2str(I)]);
            pbm.gconst(ig_(['E',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'H11'))
            pb.x0(ix_('H11'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('H11')),1) = 1.0;
        end
        if(isKey(ix_,'H12'))
            pb.x0(ix_('H12'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('H12')),1) = 0.0;
        end
        if(isKey(ix_,'H22'))
            pb.x0(ix_('H22'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('H22')),1) = 1.0;
        end
        if(isKey(ix_,'G1'))
            pb.x0(ix_('G1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('G1')),1) = 0.0;
        end
        if(isKey(ix_,'G2'))
            pb.x0(ix_('G2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('G2')),1) = 0.0;
        end
        for I=v_('1'):v_('NPTS')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_(['XD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_(['XD',int2str(I)]);
            end
            if(isKey(ix_,['Y',int2str(I)]))
                pb.x0(ix_(['Y',int2str(I)]),1) = v_(['YD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['Y',int2str(I)])),1) = v_(['YD',int2str(I)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eHXX',iet_);
        elftv{it}{1} = 'H';
        elftv{it}{2} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eHXY',iet_);
        elftv{it}{1} = 'H';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eGX',iet_);
        elftv{it}{1} = 'G';
        elftv{it}{2} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NPTS')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXX';
            ielftype(ie) = iet_('eHXX');
            vname = 'H11';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXY';
            ielftype(ie) = iet_('eHXY');
            vname = 'H12';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXX';
            ielftype(ie) = iet_('eHXX');
            vname = 'H22';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['ED',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eGX';
            ielftype(ie) = iet_('eGX');
            vname = 'G1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('G',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eGX';
            ielftype(ie) = iet_('eGX');
            vname = 'G2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('G',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NPTS')
            ig = ig_(['OX',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OY',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['E',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = 2.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['ED',int2str(I)]);
            pbm.grelw{ig}(posel) = -2.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -2.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(3)             350.29936756
% LO SOLTN(4)             1414.0524915
% LO SOLTN(5)             ???
% LO SOLTN(6)             ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QQR2-AN-V-V';
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

    case 'eHXX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(2);
            g_(2,1) = 2.0*EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = EV_(2)+EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)+EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eHXY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    case 'eGX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
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

