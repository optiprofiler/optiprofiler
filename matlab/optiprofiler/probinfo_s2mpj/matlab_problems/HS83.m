function varargout = HS83(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS83
%    *********
% 
%    Source: problem 83 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: B Baudson, Apr 1989.
% 
%    classification = 'QQR2-AN-5-3'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS83';

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
        v_('N') = 5;
        v_('A1') = 85.334407;
        v_('A2') = 0.0056858;
        v_('A3') = 0.0006262;
        v_('A4') = 0.0022053;
        v_('A5') = 80.51249;
        v_('A6') = 0.0071317;
        v_('A7') = 0.0029955;
        v_('A8') = 0.0021813;
        v_('A9') = 9.300961;
        v_('A10') = 0.0047026;
        v_('A11') = 0.0012547;
        v_('A12') = 0.0019085;
        v_('1') = 1;
        v_('CS1') = -1.0*v_('A1');
        v_('CS2') = -90.0+v_('A5');
        v_('CS2') = -1.0*v_('CS2');
        v_('CS3') = -20.0+v_('A9');
        v_('CS3') = -1.0*v_('CS3');
        v_('-A4') = -1.0*v_('A4');
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
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 37.293239+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 37.293239;
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C2';
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C3';
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
        pbm.gconst(ig_('OBJ')) = 40792.141;
        pbm.gconst(ig_('C1')) = v_('CS1');
        pbm.gconst(ig_('C2')) = v_('CS2');
        pbm.gconst(ig_('C3')) = v_('CS3');
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        grange(ig_('C1')) = 92.0;
        grange(ig_('C2')) = 20.0;
        grange(ig_('C3')) = 5.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 27.0*ones(pb.n,1);
        pb.xupper = 45.0*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 78.0;
        pb.xupper(ix_('X1')) = 102.0;
        pb.xlower(ix_('X2'),1) = 33.0;
        pb.xupper(ix_('X2')) = 45.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 27.0*ones(pb.n,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 78.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 78.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 33.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 33.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,27.0,45.0,27.0);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 5.3578547;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 0.8356891;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('-A4');
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A7');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A8');
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A10');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A11');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('A12');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -30665.53867
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QQR2-AN-5-3';
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

    case 'en2PR'

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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^2;
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

