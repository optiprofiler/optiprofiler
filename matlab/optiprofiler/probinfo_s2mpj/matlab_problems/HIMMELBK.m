function varargout = HIMMELBK(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELBK
%    *********
% 
%    A blending problem for multi-component mixtures, by Paviani.
%    It has a linear objective and linear and nonlinear constraints.
% 
%    Compared to the problem specified in Himmelblau, the inequality
%    constraints have been removed, because, as stated in this source,
%    they impose that
%    X(1)=X(2)=X(3)=X(7)=X(9)=X(9)=X(13)=X(14)=X(15)=X(19)=X(20)=X(21)=0
%    which is clearly contradictory with the given solution(s).  As
%    there does not seem to be a natural way to correct this statement
%    without knowing more about the original problem, the troublesome
%    constraints have been removed.
% 
%    Source: from problem 20 in
%    D.H. Himmelblau,
%    "Applied nonlinear programming",
%    McGraw-Hill, New-York, 1972.
% 
%    SIF input: Ph. Toint, March 1991.
% 
%    classification = 'LOR2-MN-24-14'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELBK';

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
        v_('F') = 142.22471;
        v_('B1') = 44.094;
        v_('B2') = 58.12;
        v_('B3') = 58.12;
        v_('B4') = 137.4;
        v_('B5') = 120.9;
        v_('B6') = 170.9;
        v_('B7') = 62.501;
        v_('B8') = 84.94;
        v_('B9') = 133.425;
        v_('B10') = 82.507;
        v_('B11') = 46.07;
        v_('B12') = 60.097;
        v_('B13') = 44.094;
        v_('B14') = 58.12;
        v_('B15') = 58.12;
        v_('B16') = 137.4;
        v_('B17') = 120.9;
        v_('B18') = 170.9;
        v_('B19') = 62.501;
        v_('B20') = 84.94;
        v_('B21') = 133.425;
        v_('B22') = 82.507;
        v_('B23') = 46.07;
        v_('B24') = 60.097;
        v_('C1') = 123.7;
        v_('C2') = 31.7;
        v_('C3') = 45.7;
        v_('C4') = 14.7;
        v_('C5') = 84.7;
        v_('C6') = 27.7;
        v_('C7') = 49.7;
        v_('C8') = 7.1;
        v_('C9') = 2.1;
        v_('C10') = 17.7;
        v_('C11') = 0.85;
        v_('C12') = 0.64;
        v_('D1') = 123.7;
        v_('D2') = 31.7;
        v_('D3') = 45.7;
        v_('D4') = 14.7;
        v_('D5') = 84.7;
        v_('D6') = 27.7;
        v_('D7') = 49.7;
        v_('D8') = 7.1;
        v_('D9') = 2.1;
        v_('D10') = 17.7;
        v_('D11') = 0.85;
        v_('D12') = 0.64;
        v_('1') = 1;
        v_('12') = 12;
        v_('13') = 13;
        v_('24') = 24;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for K=v_('1'):v_('24')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(K)],ix_);
            pb.xnames{iv} = ['X',int2str(K)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0693+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0693;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0577+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0577;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.05;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.2;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.26+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.26;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.55+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.55;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.06;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.1;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.12+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.12;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.18+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.18;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.1;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.09+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.09;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0693+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0693;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0577+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0577;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.05;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.2;
        end
        iv = ix_('X17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.26+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.26;
        end
        iv = ix_('X18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.55+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.55;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.06;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.1;
        end
        iv = ix_('X21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.12+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.12;
        end
        iv = ix_('X22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.18+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.18;
        end
        iv = ix_('X23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.1;
        end
        iv = ix_('X24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.09+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.09;
        end
        for I=v_('1'):v_('12')
            [ig,ig_] = s2mpjlib('ii',['CA',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CA',int2str(I)];
        end
        for I=v_('1'):v_('24')
            [ig,ig_] = s2mpjlib('ii','CA13',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CA13';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            v_('1/DI') = 1.0/v_(['D',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii','CA14',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CA14';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/DI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/DI');
            end
            v_('F/BI+12') = v_('F')/v_(['B',int2str(round(v_('I+12')))]);
            iv = ix_(['X',int2str(round(v_('I+12')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('F/BI+12')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('F/BI+12');
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
        pbm.gconst(ig_('CA13')) = 1.0;
        pbm.gconst(ig_('CA14')) = 1.671;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.04*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            for J=v_('1'):v_('12')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['X',int2str(round(v_('I+12')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for J=v_('13'):v_('24')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            for J=v_('1'):v_('12')
                v_('BI/BJ') = v_(['B',int2str(I)])/v_(['B',int2str(J)]);
                v_('40BI/BJ') = 40.0*v_('BI/BJ');
                ig = ig_(['CA',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('40BI/BJ');
            end
            for J=v_('13'):v_('24')
                v_('B+/BJ') = v_(['B',int2str(round(v_('I+12')))])/v_(['B',int2str(J)]);
                v_('CB+/BJ') = v_(['C',int2str(I)])*v_('B+/BJ');
                v_('-CB+/BJ') = -1.0*v_('CB+/BJ');
                ig = ig_(['CA',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-CB+/BJ');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                0.0893344
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MN-24-14';
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

