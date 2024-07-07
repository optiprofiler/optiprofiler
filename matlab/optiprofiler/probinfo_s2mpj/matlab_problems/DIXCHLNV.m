function varargout = DIXCHLNV(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DIXCHLNV
%    *********
% 
%    A variant of a constrained problem set as a challenge for SQP methods
%    by L.C.W. Dixon at the APMOD91 Conference.
%    The variant from the original problem is that the variables have
%    been constrained to be above 1.0D-15, which in turn allows the
%    constraints to be expressed in terms of logarithms.
% 
%    Source: a modification (by Ph. Toint) of
%    L.C.W. Dixon, personnal communication, Jan 1991.
% 
%    SIF input: Ph. Toint, Feb 1991.
% 
%    classification = 'SOR2-AN-V-V'
% 
%    Number of variables
%    (variable, but must be even and at least equal to 4)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER     original value
% IE N                   1000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DIXCHLNV';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5000           $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('N-3') = -3+v_('N');
        v_('N-1') = -1+v_('N');
        v_('90.0') = 90.0;
        v_('10.1') = 10.1;
        v_('19.8') = 19.8;
        v_('1/90.0') = 1.0/v_('90.0');
        v_('1/10.1') = 1.0/v_('10.1');
        v_('1/19.8') = 1.0/v_('19.8');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-3')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            [ig,ig_] = s2mpjlib('ii',['A',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = 0.01;
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = v_('1/90.0');
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = v_('1/10.1');
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = v_('1/10.1');
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('1/19.8');
        end
        for I=v_('2'):v_('2'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['P',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['P',int2str(I)];
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
        for I=v_('1'):v_('N-3')
            pbm.gconst(ig_(['B',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['D',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['E',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['F',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 1.0e-8*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('X0') = 2.0;
        for I=v_('1'):v_('2'):v_('N-1')
            pb.x0(ix_(['X',int2str(I)]),1) = v_('X0');
            v_('1/X0') = 1.0/v_('X0');
            v_('I+1') = 1+I;
            pb.x0(ix_(['X',int2str(round(v_('I+1')))]),1) = v_('1/X0');
            v_('X0') = 1.0+v_('X0');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eS2PR',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eLOGV',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N-1')
            ename = ['XSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('N-3')
            v_('I+1') = 1+I;
            v_('I+3') = 3+I;
            ename = ['PR',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eS2PR';
            ielftype(ie) = iet_('eS2PR');
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('N')
            ename = ['LOGX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eLOGV';
            ielftype(ie) = iet_('eLOGV');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,1.0e-8,[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N-3')
            v_('I+2') = 2+I;
            ig = ig_(['A',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['B',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['C',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(round(v_('I+2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['D',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['E',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['F',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PR',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for I=v_('2'):v_('2'):v_('N')
            for J=v_('1'):I
                ig = ig_(['P',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['LOGX',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
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
        pb.pbclass = 'SOR2-AN-V-V';
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eLOGV'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = log(EV_(1));
        if(nargout>1)
            g_(1,1) = 1.0/EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -1.0/(EV_(1)*EV_(1));
                varargout{3} = H_;
            end
        end

    case 'eS2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-1.0)*(EV_(2)-1.0);
        if(nargout>1)
            g_(1,1) = EV_(2)-1.0;
            g_(2,1) = EV_(1)-1.0;
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

