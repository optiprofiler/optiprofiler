function varargout = HS56(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS56
%    *********
% 
%    Source: problem 56 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: A.R. Conn, April 1990
% 
%    classification = 'OOR2-AN-7-4'
% 
%    some useful parameters, including N, the number of variables.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS56';

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
        v_('N') = 7;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_(['PAR',int2str(round(v_('1')))]) = 4.2;
        v_(['PAR',int2str(round(v_('2')))]) = 4.2;
        v_(['PAR',int2str(round(v_('3')))]) = 4.2;
        v_(['PAR',int2str(round(v_('4')))]) = 7.2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('1')))];
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('2')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('2')))];
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('3')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('3')))];
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('4')))];
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        [ig,ig_] = s2mpjlib('ii',['CON',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['CON',int2str(round(v_('4')))];
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
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
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 0.50973968;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 0.50973968;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 0.50973968;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 0.50973968;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 0.50973968;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 0.50973968;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 0.98511078;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 0.98511078;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en3PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        [it,iet_] = s2mpjlib( 'ii', 'ePSNSQ',iet_);
        elftv{it}{1} = 'V1';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PROD';
        ielftype(ie) = iet_('en3PROD');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('1'):v_('4')
            v_('J') = 1+I;
            v_('K') = 3+I;
            ename = ['E',int2str(round(v_('J')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePSNSQ';
            ielftype(ie) = iet_('ePSNSQ');
            ename = ['E',int2str(round(v_('J')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('K')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('J')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['PAR',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        for I=v_('1'):v_('4')
            v_('J') = 1+I;
            ig = ig_(['CON',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('J')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                 -3.456
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-7-4';
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

    case 'en3PROD'

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

    case 'ePSNSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SUM = EV_(1)+EV_(1);
        varargout{1} = -pbm.elpar{iel_}(1)*sin(EV_(1))^2;
        if(nargout>1)
            g_(1,1) = -pbm.elpar{iel_}(1)*sin(SUM);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -2.0*pbm.elpar{iel_}(1)*cos(SUM);
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

