function varargout = PFIT4LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PFIT4LS
%    *********
% 
%    The problem is to fit a model containing a pole, given data
%    for values, first and second derivatives at two distinct points.
%    This is a least-squares version of problem PFIT4.
% 
%    The problem is not convex.
% 
%    SIF input: Ph. Toint, March 1994.
%               Lower bound on H added, Nov 2002.
% 
%    classification = 'SBR2-AN-3-0'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PFIT4LS';

switch(action)

    case 'setup'

    pb.name      = 'PFIT4LS';
    pb.sifpbname = 'PFIT4LS';
    pbm.name     = 'PFIT4LS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('CF') = -98.96296296;
        v_('CG') = -216.0987654;
        v_('CH') = -239.6707818;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','A',ix_);
        pb.xnames{iv} = 'A';
        [iv,ix_] = s2xlib('ii','R',ix_);
        pb.xnames{iv} = 'R';
        [iv,ix_] = s2xlib('ii','H',ix_);
        pb.xnames{iv} = 'H';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','EF',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','EG',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','EH',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('EF')) = v_('CF');
        pbm.gconst(ig_('EG')) = v_('CG');
        pbm.gconst(ig_('EH')) = v_('CH');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_('H'),1) = -0.5;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('A'),1) = 1.0;
        pb.x0(ix_('R'),1) = 0.0;
        pb.x0(ix_('H'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'T1',iet_);
        elftv{it}{1} = 'AA';
        elftv{it}{2} = 'RR';
        elftv{it}{3} = 'XX';
        [it,iet_] = s2xlib( 'ii', 'T2',iet_);
        elftv{it}{1} = 'AA';
        elftv{it}{2} = 'RR';
        elftv{it}{3} = 'XX';
        [it,iet_] = s2xlib( 'ii', 'T3',iet_);
        elftv{it}{1} = 'AA';
        elftv{it}{2} = 'RR';
        elftv{it}{3} = 'XX';
        [it,iet_] = s2xlib( 'ii', 'T4',iet_);
        elftv{it}{1} = 'AA';
        elftv{it}{2} = 'RR';
        elftv{it}{3} = 'XX';
        [it,iet_] = s2xlib( 'ii', 'T5',iet_);
        elftv{it}{1} = 'AA';
        elftv{it}{2} = 'RR';
        elftv{it}{3} = 'XX';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'EA';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'T3';
        ielftype(ie) = iet_('T3');
        vname = 'A';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RR',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'H';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'T2';
        ielftype(ie) = iet_('T2');
        vname = 'A';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RR',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'H';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EC';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'T1';
        ielftype(ie) = iet_('T1');
        vname = 'A';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RR',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'H';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'ED';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'T4';
        ielftype(ie) = iet_('T4');
        vname = 'A';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RR',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'H';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EE';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'T5';
        ielftype(ie) = iet_('T5');
        vname = 'A';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RR',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'H';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('XX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','L2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'L2';
        end
        ig = ig_('EF');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EA');
        pbm.grelw{ig}(posel) = -0.5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EC');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ED');
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('EG');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EA');
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EB');
        pbm.grelw{ig}(posel) = 1;
        ig = ig_('EH');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EE');
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = 'SBR2-AN-3-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'T1'

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

    case 'T2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A1 = EV_(1)+1.0;
        Y = 1.0+EV_(3);
        LOGY = log(Y);
        C = Y^(-A1);
        CC = C/Y;
        CCC = CC/Y;
        B = 1.0-C;
        BA = LOGY*C;
        BX = A1*CC;
        BAA = -LOGY*LOGY*C;
        BAX = -LOGY*BX+CC;
        BXX = -A1*(A1+1.0)*CCC;
        ARX = EV_(1)*EV_(2)*EV_(3);
        varargout{1} = ARX*B;
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*B+ARX*BA;
            g_(2,1) = EV_(1)*EV_(3)*B;
            g_(3,1) = EV_(1)*EV_(2)*B+ARX*BX;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0*EV_(2)*EV_(3)*BA+ARX*BAA;
                H_(1,2) = EV_(3)*B+EV_(1)*EV_(3)*BA;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*B+EV_(2)*EV_(3)*BX+EV_(1)*EV_(2)*BA+ARX*BAX;
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1)*B+EV_(1)*EV_(3)*BX;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*EV_(1)*EV_(2)*BX+ARX*BXX;
                varargout{3} = H_;
            end
        end

    case 'T3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*(EV_(1)+1.0)*EV_(2)*EV_(3)*EV_(3);
        if(nargout>1)
            g_(1,1) = (2.0*EV_(1)+1.0)*EV_(2)*EV_(3)*EV_(3);
            g_(2,1) = EV_(1)*(EV_(1)+1.0)*EV_(3)*EV_(3);
            g_(3,1) = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(2)*EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0*EV_(2)*EV_(3)*EV_(3);
                H_(1,2) = (2.0*EV_(1)+1.0)*EV_(3)*EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0*(2.0*EV_(1)+1.0)*EV_(2)*EV_(3);
                H_(3,1) = H_(1,3);
                H_(2,3) = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(3);
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(2);
                varargout{3} = H_;
            end
        end

    case 'T4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        Y = 1.0+EV_(3);
        LOGY = log(Y);
        C = Y^(-EV_(1));
        CC = C/Y;
        CCC = CC/Y;
        B = 1.0-C;
        BA = LOGY*C;
        BX = EV_(1)*CC;
        BAA = -LOGY*LOGY*C;
        BAX = -LOGY*BX+CC;
        BXX = -EV_(1)*(EV_(1)+1.0)*CCC;
        varargout{1} = EV_(2)*B;
        if(nargout>1)
            g_(1,1) = EV_(2)*BA;
            g_(2,1) = B;
            g_(3,1) = EV_(2)*BX;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = EV_(2)*BAA;
                H_(1,2) = BA;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*BAX;
                H_(3,1) = H_(1,3);
                H_(2,3) = BX;
                H_(3,2) = H_(2,3);
                H_(3,3) = EV_(2)*BXX;
                varargout{3} = H_;
            end
        end

    case 'T5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A1 = EV_(1)+2.0;
        Y = 1.0+EV_(3);
        LOGY = log(Y);
        C = Y^(-A1);
        CC = C/Y;
        CCC = CC/Y;
        B = 1.0-C;
        BA = LOGY*C;
        BX = A1*CC;
        BAA = -LOGY*LOGY*C;
        BAX = -LOGY*BX+CC;
        BXX = -A1*(A1+1.0)*CCC;
        D = EV_(1)*(EV_(1)+1.0)*EV_(2)*EV_(3)*EV_(3);
        DA = (2.0*EV_(1)+1.0)*EV_(2)*EV_(3)*EV_(3);
        DR = EV_(1)*(EV_(1)+1.0)*EV_(3)*EV_(3);
        DX = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(2)*EV_(3);
        DAA = 2.0*EV_(2)*EV_(3)*EV_(3);
        DAR = (2.0*EV_(1)+1.0)*EV_(3)*EV_(3);
        DAX = 2.0*(2.0*EV_(1)+1.0)*EV_(2)*EV_(3);
        DRX = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(3);
        DXX = 2.0*EV_(1)*(EV_(1)+1.0)*EV_(2);
        varargout{1} = D*B;
        if(nargout>1)
            g_(1,1) = DA*B+D*BA;
            g_(2,1) = DR*B;
            g_(3,1) = DX*B+D*BX;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = DAA*B+2.0*DA*BA+D*BAA;
                H_(1,2) = DAR*B+DR*BA;
                H_(2,1) = H_(1,2);
                H_(1,3) = DAX*B+DA*BX+DX*BA+D*BAX;
                H_(3,1) = H_(1,3);
                H_(2,3) = DRX*B+DR*BX;
                H_(3,2) = H_(2,3);
                H_(3,3) = DXX*B+2.0*DX*BX+D*BXX;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'L2'

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
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
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

