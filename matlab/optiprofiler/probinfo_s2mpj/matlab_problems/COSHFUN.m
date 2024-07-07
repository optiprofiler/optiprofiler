function varargout = COSHFUN(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : COSHFUN
%    *********
% 
%    A nonlinear minmax problem.
% 
%    Source:
%    K. Jonasson and K. Madsen,
%    "Corrected sequential linear programming for sparse
%    minimax optimization", Technical report, Institute for Numerical
%    Analysis, Technical U. of Denmark.
% 
%    SIF input: Nick Gould, October 1992.
% 
%    classification = 'LOR2-AN-V-V'
% 
%   the number of functions
% 
%       Alternative values for the SIF file parameters:
% IE M                   3              $-PARAMETER
% IE M                   8              $-PARAMETER
% IE M                   14             $-PARAMETER
% IE M                   20             $-PARAMETER     original value
% IE M                   200            $-PARAMETER
% IE M                   2000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'COSHFUN';

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
            v_('M') = 8;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
        v_('N') = 3*v_('M');
        v_('N-3') = -3+v_('N');
        v_('N-5') = -5+v_('N');
        v_('0') = 0;
        v_('1') = 1;
        v_('3') = 3;
        v_('6') = 6;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','F',ix_);
        pb.xnames{iv} = 'F';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('F');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(I)];
            iv = ix_('F');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C1';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('6'):v_('3'):v_('N-3')
            v_('I-5') = -5+I;
            v_('I+3') = 3+I;
            v_('I/3') = fix(I/v_('3'));
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('I/3')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('I/3')))];
            iv = ix_(['X',int2str(round(v_('I-5')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0;
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('I/3')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('I/3')))];
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('M')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('M')))];
        iv = ix_(['X',int2str(round(v_('N-5')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eCOSH',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('3'):v_('3'):v_('N')
            v_('I-1') = -1+I;
            v_('I-2') = -2+I;
            v_('I/3') = fix(I/v_('3'));
            ename = ['SQR',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            ename = ['SQR',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['COSH',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOSH';
            ielftype(ie) = iet_('eCOSH');
            ename = ['COSH',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['PROD',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            ename = ['PROD',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I-2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['PROD',int2str(round(v_('I/3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['C',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['SQR',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['COSH',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PROD',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0e+0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0;
                varargout{3} = H_;
            end
        end

    case 'eCOSH'

        EV_  = varargin{1};
        iel_ = varargin{2};
        COSHX = cosh(EV_(1));
        varargout{1} = COSHX;
        if(nargout>1)
            g_(1,1) = sinh(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = COSHX;
                varargout{3} = H_;
            end
        end

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 2.0e+0*EV_(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 4.0e+0*EV_(1)*EV_(2);
            g_(2,1) = 2.0e+0*EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 4.0e+0*EV_(2);
                H_(1,2) = 4.0e+0*EV_(1);
                H_(2,1) = H_(1,2);
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

