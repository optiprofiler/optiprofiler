function varargout = CSFI2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CSFI2
%    *********
% 
%    Source: problem MINLEN in
%    Vasko and Stott
%    "Optimizing continuous caster product dimensions:
%     an example of a nonlinear design problem in the steel industry"
%    SIAM Review, Vol 37 No, 1 pp.82-84, 1995
% 
%    SIF input: A.R. Conn April 1995
% 
%    classification = 'LOR2-RN-5-4'
% 
%    input parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CSFI2';

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
        v_('MINTPH') = 45.0;
        v_('MINTHICK') = 7.0;
        v_('MINAREA') = 200.0;
        v_('MAXAREA') = 250.0;
        v_('MAXASPR') = 2.0;
        v_('K') = 1.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','THICK',ix_);
        pb.xnames{iv} = 'THICK';
        [iv,ix_] = s2mpjlib('ii','WID',ix_);
        pb.xnames{iv} = 'WID';
        [iv,ix_] = s2mpjlib('ii','LEN',ix_);
        pb.xnames{iv} = 'LEN';
        [iv,ix_] = s2mpjlib('ii','TPH',ix_);
        pb.xnames{iv} = 'TPH';
        [iv,ix_] = s2mpjlib('ii','IPM',ix_);
        pb.xnames{iv} = 'IPM';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('LEN');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CIPM',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CIPM';
        iv = ix_('IPM');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CLEN',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CLEN';
        iv = ix_('LEN');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','WOT',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'WOT';
        [ig,ig_] = s2mpjlib('ii','TTW',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'TTW';
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
        pbm.gconst(ig_('WOT')) = v_('MAXASPR');
        pbm.gconst(ig_('TTW')) = v_('MINAREA');
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        v_('RHS') = v_('MAXAREA')-v_('MINAREA');
        grange(ig_('TTW')) = v_('RHS');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('THICK'),1) = v_('MINTHICK');
        pb.xlower(ix_('TPH'),1) = v_('MINTPH');
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eCMPLQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        [it,iet_] = s2mpjlib( 'ii', 'eSQQUT',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'eQUOTN',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCMPLQ';
        ielftype(ie) = iet_('eCMPLQ');
        vname = 'TPH';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'WID';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'THICK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQUT';
        ielftype(ie) = iet_('eSQQUT');
        vname = 'THICK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'IPM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQUOTN';
        ielftype(ie) = iet_('eQUOTN');
        vname = 'WID';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'THICK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'THICK';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'WID';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('CIPM');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('CLEN');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('WOT');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('TTW');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               55.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-RN-5-4';
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

    case 'eCMPLQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP0 = EV_(2)*EV_(3);
        TMP1 = 117.3708920187793427e0*EV_(1)/TMP0;
        TMP2 = 117.3708920187793427e0/TMP0;
        varargout{1} = TMP1;
        if(nargout>1)
            g_(1,1) = TMP2;
            g_(2,1) = -TMP1/EV_(2);
            g_(3,1) = -TMP1/EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = -TMP2/EV_(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = -TMP2/EV_(3);
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0e0*TMP1/(EV_(2)*EV_(2));
                H_(2,3) = TMP1/TMP0;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0e0*TMP1/(EV_(3)*EV_(3));
                varargout{3} = H_;
            end
        end

    case 'eSQQUT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP = EV_(1)*EV_(2)/48.0e0;
        varargout{1} = EV_(1)*TMP;
        if(nargout>1)
            g_(1,1) = 2.0e0*TMP;
            g_(2,1) = EV_(1)*EV_(1)/48.0e0;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = EV_(2)/24.0e0;
                H_(1,2) = EV_(1)/24.0e0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eQUOTN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMP = EV_(1)/EV_(2);
        varargout{1} = TMP;
        if(nargout>1)
            g_(1,1) = 1.0e0/EV_(2);
            g_(2,1) = -TMP/EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0e0/(EV_(2)*EV_(2));
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0e0*TMP/(EV_(2)*EV_(2));
                varargout{3} = H_;
            end
        end

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e0;
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

