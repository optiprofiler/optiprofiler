function varargout = ALLINITC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ALLINITC
%    *********
% 
%    A problem with "all in it". Intended to verify that changes
%    to LANCELOT are safe.
% 
%    Source:
%    N. Gould: private communication.
% 
%    SIF input: Nick Gould, June 1990.
% 
%    classification = 'C-COOR2-AY-4-1'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ALLINITC';

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
        [ig,ig_] = s2mpjlib('ii','FT1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FT2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FT3',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FT4',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FT5',ig_);
        gtype{ig} = '<>';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FT6',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FNT1',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FNT2',ig_);
        gtype{ig} = '<>';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FNT3',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FNT4',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','FNT5',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FNT6',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
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
        pbm.gconst(ig_('FT2')) = 1.0;
        pbm.gconst(ig_('FT5')) = 3.0;
        pbm.gconst(ig_('FNT2')) = 1.0;
        pbm.gconst(ig_('FNT5')) = 4.0;
        pbm.gconst(ig_('C1')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1')) = -Inf;
        pb.xupper(ix_('X1'),1) = +Inf;
        pb.xlower(ix_('X2'),1) = 1.0;
        pb.xlower(ix_('X3'),1) = -1.0e+10;
        pb.xupper(ix_('X3')) = 1.0;
        pb.xlower(ix_('X4'),1) = 2.0;
        pb.xupper(ix_('X4'),1) = 2.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR2',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSINSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'ePRODSQR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'FT3E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT4E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT4E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR2';
        ielftype(ie) = iet_('eSQR2');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT56E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINSQR';
        ielftype(ie) = iet_('eSINSQR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FT5E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODSQR';
        ielftype(ie) = iet_('ePRODSQR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT3E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT4E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT4E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR2';
        ielftype(ie) = iet_('eSQR2');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT56E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSINSQR';
        ielftype(ie) = iet_('eSINSQR');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FNT5E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODSQR';
        ielftype(ie) = iet_('ePRODSQR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gTRIVIAL',igt_);
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_('FT1');
        pbm.grftype{ig} = 'gTRIVIAL';
        ig = ig_('FT2');
        pbm.grftype{ig} = 'gTRIVIAL';
        ig = ig_('FT3');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT3E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT4');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT4E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FT4E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT5');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT56E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FT5E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FT6');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT56E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT3E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT4E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FNT4E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT56E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FNT5E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('FNT6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FNT56E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('C1');
        pbm.grftype{ig} = 'gTRIVIAL';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FT3E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('FT4E1');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AY-4-1';
        pb.x0          = zeros(pb.n,1);
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

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

    case 'eSQR2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eSINSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SINX = sin(EV_(1));
        COSX = cos(EV_(1));
        varargout{1} = SINX*SINX;
        if(nargout>1)
            g_(1,1) = 2.0*SINX*COSX;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*(COSX*COSX-SINX*SINX);
                varargout{3} = H_;
            end
        end

    case 'ePRODSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XX = EV_(1)*EV_(1);
        YY = EV_(2)*EV_(2);
        varargout{1} = XX*YY;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*YY;
            g_(2,1) = 2.0*XX*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*YY;
                H_(1,2) = 4.0*EV_(1)*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*XX;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gTRIVIAL'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_;
        if(nargout>1)
            g_ = 1.0;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 0.0;
                varargout{3} = H_;
            end
        end

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

