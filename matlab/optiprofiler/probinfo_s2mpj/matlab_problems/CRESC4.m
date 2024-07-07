function varargout = CRESC4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CRESC4
%    *********
% 
%    This problem consists in finding the crescent of smallest area containing
%    a set of points given in the plane.   This problem arises as a subproblem
%    in pattern recognition and has been suggested by J.P. Rasson.  It
%    originates in the detection of "salt domes" (with the potential presence of
%    oil pockets!) from geological data.
% 
%    The present problem is a simplified version where the crescent is entirely
%    determined by the only four data points.
% 
%    The problem is not convex.
% 
%    A crescent is defined as follows.  Assume one has two circles of center
%    C1 and C2 and of radii r1 and r2 respectively. Assume furthermore that
%    r1 >= r2 and that C2 is within C1.  Assume finally that the distance from
%    C1 to C2 is >= r1 - r2.  Then the crescent is the part of the plane
%    contained in circle 2 but not in circle 1.
% 
%    In order to preserve feasibility at all stages (ensuring that the
%    crescent exists and that its area can be computed), the following
%    parametrization is used:
% 
%    ( C2x, C2y ) = ( C1x, C1y ) + a * d * ( cos(t), sin(t) )
% 
%    r1 = a * d + r 
% 
%    r2 = ( a + 1 ) * d + r
% 
%    with the bounds
% 
%    a >= 1, 0 <= t <= 2 * pi, r2 >= 0 , 0 <= d <= 1.
% 
%    SIF input: Ph. Toint, June 1993.
% 
%    classification = 'OOR2-MY-6-8'
% 
%    number of points to be included in the crescent.
%    the number of constraints is 2*NP
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CRESC4';

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
        v_('NP') = 4;
        v_('X1') = 1.0;
        v_('Y1') = 0.0;
        v_('X2') = 0.0;
        v_('Y2') = 1.0;
        v_('X3') = 0.0;
        v_('Y3') = -1.0;
        v_('X4') = 0.5;
        v_('Y4') = 0.0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','V1',ix_);
        pb.xnames{iv} = 'V1';
        [iv,ix_] = s2mpjlib('ii','W1',ix_);
        pb.xnames{iv} = 'W1';
        [iv,ix_] = s2mpjlib('ii','D',ix_);
        pb.xnames{iv} = 'D';
        [iv,ix_] = s2mpjlib('ii','A',ix_);
        pb.xnames{iv} = 'A';
        [iv,ix_] = s2mpjlib('ii','T',ix_);
        pb.xnames{iv} = 'T';
        [iv,ix_] = s2mpjlib('ii','R',ix_);
        pb.xnames{iv} = 'R';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('NP')
            [ig,ig_] = s2mpjlib('ii',['IS2',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['IS2',int2str(I)];
            [ig,ig_] = s2mpjlib('ii',['OS1',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['OS1',int2str(I)];
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('V1')) = -Inf;
        pb.xupper(ix_('V1'),1) = +Inf;
        pb.xlower(ix_('W1')) = -Inf;
        pb.xupper(ix_('W1'),1) = +Inf;
        pb.xlower(ix_('R'),1) = 0.39;
        pb.xlower(ix_('A'),1) = 1.0;
        pb.xlower(ix_('T'),1) = 0.0;
        pb.xupper(ix_('T')) = 6.2831852;
        pb.xlower(ix_('D'),1) = 1.0e-8;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('V1'),1) = -40.0;
        pb.x0(ix_('W1'),1) = 5.0;
        pb.x0(ix_('R'),1) = 0.75;
        pb.x0(ix_('A'),1) = 2.0;
        pb.x0(ix_('T'),1) = 1.5;
        pb.x0(ix_('D'),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR1',iet_);
        elftv{it}{1} = 'A';
        elftv{it}{2} = 'R';
        elftv{it}{3} = 'D';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR2',iet_);
        elftv{it}{1} = 'D';
        elftv{it}{2} = 'R';
        [it,iet_] = s2mpjlib( 'ii', 'eSC',iet_);
        elftv{it}{1} = 'AZ';
        elftv{it}{2} = 'BZ';
        elftv{it}{3} = 'DZ';
        [it,iet_] = s2mpjlib( 'ii', 'eDIST',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eDISTX',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'T';
        elftv{it}{3} = 'A';
        elftv{it}{4} = 'D';
        elftp{it}{1} = 'P';
        [it,iet_] = s2mpjlib( 'ii', 'eDISTY',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'T';
        elftv{it}{3} = 'A';
        elftv{it}{4} = 'D';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OB';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSC';
        ielftype(ie) = iet_('eSC');
        vname = 'A';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('AZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('BZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('DZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'R2SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR2';
        ielftype(ie) = iet_('eSQR2');
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('D',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('R',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'R1SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR1';
        ielftype(ie) = iet_('eSQR1');
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('D',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'A';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('A',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'R';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('R',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('1'):v_('NP')
            ename = ['XV1',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDIST';
            ielftype(ie) = iet_('eDIST');
            vname = 'V1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
            ename = ['XV2',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDISTX';
            ielftype(ie) = iet_('eDISTX');
            vname = 'V1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'A';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'D';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('D',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'T';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
            ename = ['YW1',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDIST';
            ielftype(ie) = iet_('eDIST');
            vname = 'W1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
            ename = ['YW2',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDISTY';
            ielftype(ie) = iet_('eDISTY');
            vname = 'W1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'A';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'D';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('D',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'T';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('P',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OB');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('1'):v_('NP')
            ig = ig_(['IS2',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XV2',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['YW2',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_('R2SQ');
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            ig = ig_(['OS1',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XV1',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['YW1',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_('R1SQ');
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution             0.87189692
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MY-6-8';
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

    case 'eSQR1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        Q = EV_(1)*EV_(3)+EV_(2);
        varargout{1} = Q*Q;
        if(nargout>1)
            g_(1,1) = 2.0*Q*EV_(3);
            g_(3,1) = 2.0*Q*EV_(1);
            g_(2,1) = 2.0*Q;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 2.0*EV_(3)*EV_(3);
                H_(1,3) = 2.0*(EV_(1)*EV_(3)+Q);
                H_(3,1) = H_(1,3);
                H_(1,2) = 2.0*EV_(3);
                H_(2,1) = H_(1,2);
                H_(3,3) = 2.0*EV_(1)*EV_(1);
                H_(3,2) = 2.0*EV_(1);
                H_(2,3) = H_(3,2);
                H_(2,2) = 2.0;
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

    case 'eSC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = EV_(3)+EV_(2);
        AB = 1.0;
        AD = 1.0;
        B = EV_(1)*EV_(3)+EV_(2);
        BA = EV_(3);
        BB = 1.0;
        BD = EV_(1);
        BAD = 1.0;
        D = EV_(1)*EV_(3);
        DA = EV_(3);
        DD = EV_(1);
        DAD = 1.0;
        E = 2.0*A*D;
        EA = 2.0*A*DA;
        EB = 2.0*AB*D;
        ED = 2.0*(A*DD+AD*D);
        EAB = 2.0*AB*DA;
        EAD = 2.0*(AD*DA+A*DAD);
        EBD = 2.0*AB*DD;
        EDD = 2.0*(AD*DD+AD*DD);
        P = 2.0*B*D;
        PA = 2.0*(B*DA+BA*D);
        loc_PB = 2.0*BB*D;
        PD = 2.0*(B*DD+BD*D);
        PAA = 2.0*(BA*DA+BA*DA);
        PAB = 2.0*BB*DA;
        PAD = 2.0*(BD*DA+B*DAD+BAD*D+BA*DD);
        PBD = 2.0*BB*DD;
        PDD = 2.0*(BD*DD+BD*DD);
        F = D*D+B*B-A*A;
        FA = 2.0*(D*DA+B*BA);
        FB = 2.0*(B*BB-A*AB);
        FD = 2.0*(D*DD+B*BD-A*AD);
        FAA = 2.0*(DA*DA+BA*BA);
        FAB = 2.0*BB*BA;
        FAD = 2.0*(DD*DA+D*DAD+BD*BA+B*BAD);
        FBB = 2.0*(BB*BB-AB*AB);
        FBD = 2.0*(BD*BB-AD*AB);
        FDD = 2.0*(DD*DD+BD*BD-AD*AD);
        G = D*D-B*B+A*A;
        GA = 2.0*(D*DA-B*BA);
        GB = 2.0*(-B*BB+A*AB);
        GD = 2.0*(D*DD-B*BD+A*AD);
        GAA = 2.0*(DA*DA-BA*BA);
        GAB = -2.0*BB*BA;
        GAD = 2.0*(DD*DA+D*DAD-BD*BA-B*BAD);
        GBB = 2.0*(-BB*BB+AB*AB);
        GBD = 2.0*(-BD*BB+AD*AB);
        GDD = 2.0*(DD*DD-BD*BD+AD*AD);
        H = F/P;
        I = FA-H*PA;
        J = FB-H*loc_PB;
        K = FD-H*PD;
        HA = I/P;
        HB = J/P;
        HD = K/P;
        IA = FAA-HA*PA-H*PAA;
        IB = FAB-HB*PA-H*PAB;
        ID = FAD-HD*PA-H*PAD;
        JB = FBB-HB*loc_PB;
        JD = FBD-HD*loc_PB-H*PBD;
        KD = FDD-HD*PD-H*PDD;
        HAA = (IA-HA*PA)/P;
        HAB = (IB-HA*loc_PB)/P;
        HAD = (ID-HA*PD)/P;
        HBB = (JB-HB*loc_PB)/P;
        HBD = (JD-HB*PD)/P;
        HDD = (KD-HD*PD)/P;
        L = -G/E;
        M = -GA-L*EA;
        N = -GB-L*EB;
        O = -GD-L*ED;
        LA = M/E;
        LB = N/E;
        LD = O/E;
        MA = -GAA-LA*EA;
        MB = -GAB-LB*EA-L*EAB;
        MD = -GAD-LD*EA-L*EAD;
        NB = -GBB-LB*EB;
        ND = -GBD-LD*EB-L*EBD;
        OD = -GDD-LD*ED-L*EDD;
        LAA = (MA-LA*EA)/E;
        LAB = (MB-LA*EB)/E;
        LAD = (MD-LA*ED)/E;
        LBB = (NB-LB*EB)/E;
        LBD = (ND-LB*ED)/E;
        LDD = (OD-LD*ED)/E;
        C = acos(H);
        CH = -1.0/sqrt(1.0-H*H);
        CHH = CH*H/(1.0-H*H);
        CA = CH*HA;
        CB = CH*HB;
        CD = CH*HD;
        CAA = CHH*HA*HA+CH*HAA;
        CAB = CHH*HA*HB+CH*HAB;
        CAD = CHH*HA*HD+CH*HAD;
        CBB = CHH*HB*HB+CH*HBB;
        CBD = CHH*HB*HD+CH*HBD;
        CDD = CHH*HD*HD+CH*HDD;
        Q = acos(L);
        QL = -1.0/sqrt(1.0-L*L);
        QLL = QL*L/(1.0-L*L);
        QA = QL*LA;
        QB = QL*LB;
        QD = QL*LD;
        QAA = QLL*LA*LA+QL*LAA;
        QAB = QLL*LA*LB+QL*LAB;
        QAD = QLL*LA*LD+QL*LAD;
        QBB = QLL*LB*LB+QL*LBB;
        QBD = QLL*LB*LD+QL*LBD;
        QDD = QLL*LD*LD+QL*LDD;
        R = B*B*C;
        RA = 2.0*B*BA*C+B*B*CA;
        RB = 2.0*B*BB*C+B*B*CB;
        RD = 2.0*B*BD*C+B*B*CD;
        RAA = 2.0*(BA*BA*C+B*BA*CA+B*BA*CA)+B*B*CAA;
        RAB = 2.0*(BB*BA*C+B*BA*CB+B*BB*CA)+B*B*CAB;
        RAD = 2.0*(BD*BA*C+B*BAD*C+B*BA*CD+B*BD*CA)+B*B*CAD;
        RBB = 2.0*(BB*BB*C+B*BB*CB+B*BB*CB)+B*B*CBB;
        RBD = 2.0*(BD*BB*C+B*BB*CD+B*BD*CB)+B*B*CBD;
        RDD = 2.0*(BD*BD*C+B*BD*CD+B*BD*CD)+B*B*CDD;
        S = A*A*Q;
        SA = A*A*QA;
        SB = 2.0*A*AB*Q+A*A*QB;
        SD = 2.0*A*AD*Q+A*A*QD;
        SAA = A*A*QAA;
        SAB = 2.0*A*AB*QA+A*A*QAB;
        SAD = 2.0*A*AD*QA+A*A*QAD;
        SBB = 2.0*(AB*AB*Q+A*AB*QB+A*AB*QB)+A*A*QBB;
        SBD = 2.0*(AD*AB*Q+A*AB*QD+A*AD*QB)+A*A*QBD;
        SDD = 2.0*(AD*AD*Q+A*AD*QD+A*AD*QD)+A*A*QDD;
        SQ = sin(Q);
        CQ = L;
        W = 0.5*E*SQ;
        WA = 0.5*(EA*SQ+E*CQ*QA);
        WB = 0.5*(EB*SQ+E*CQ*QB);
        WD = 0.5*(ED*SQ+E*CQ*QD);
        WAA = 0.5*(EA*CQ*QA+EA*CQ*QA-E*SQ*QA*QA+E*CQ*QAA);
        WAB = 0.5*(EAB*SQ+EA*CQ*QB+EB*CQ*QA-E*SQ*QB*QA+E*CQ*QAB);
        WAD = 0.5*(EAD*SQ+EA*CQ*QD+ED*CQ*QA-E*SQ*QD*QA+E*CQ*QAD);
        WBB = 0.5*(EB*CQ*QB+EB*CQ*QB-E*SQ*QB*QB+E*CQ*QBB);
        WBD = 0.5*(EBD*SQ+EB*CQ*QD+ED*CQ*QB-E*SQ*QD*QB+E*CQ*QBD);
        WDD = 0.5*(EDD*SQ+ED*CQ*QD+ED*CQ*QD-E*SQ*QD*QD+E*CQ*QDD);
        V = S-R+W;
        VA = SA-RA+WA;
        VB = SB-RB+WB;
        VD = SD-RD+WD;
        VAA = SAA-RAA+WAA;
        VAB = SAB-RAB+WAB;
        VAD = SAD-RAD+WAD;
        VBB = SBB-RBB+WBB;
        VBD = SBD-RBD+WBD;
        VDD = SDD-RDD+WDD;
        varargout{1} = V;
        if(nargout>1)
            g_(1,1) = VA;
            g_(2,1) = VB;
            g_(3,1) = VD;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = VAA;
                H_(1,2) = VAB;
                H_(2,1) = H_(1,2);
                H_(1,3) = VAD;
                H_(3,1) = H_(1,3);
                H_(2,2) = VBB;
                H_(2,3) = VBD;
                H_(3,2) = H_(2,3);
                H_(3,3) = VDD;
                varargout{3} = H_;
            end
        end

    case 'eDIST'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-pbm.elpar{iel_}(1))^2;
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-pbm.elpar{iel_}(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eDISTY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ST = sin(EV_(2));
        CT = cos(EV_(2));
        B = EV_(1)+EV_(3)*EV_(4)*ST-pbm.elpar{iel_}(1);
        BA = EV_(4)*ST;
        BD = EV_(3)*ST;
        BT = EV_(3)*EV_(4)*CT;
        varargout{1} = B*B;
        if(nargout>1)
            g_(1,1) = 2.0*B;
            g_(3,1) = 2.0*B*BA;
            g_(4,1) = 2.0*B*BD;
            g_(2,1) = 2.0*B*BT;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 2.0;
                H_(1,3) = 2.0*BA;
                H_(3,1) = H_(1,3);
                H_(1,4) = 2.0*BD;
                H_(4,1) = H_(1,4);
                H_(1,2) = 2.0*BT;
                H_(2,1) = H_(1,2);
                H_(3,3) = 2.0*BA*BA;
                H_(3,4) = 2.0*(BD*BA+B*ST);
                H_(4,3) = H_(3,4);
                H_(3,2) = 2.0*(BT*BA+B*EV_(4)*CT);
                H_(2,3) = H_(3,2);
                H_(4,4) = 2.0*BD*BD;
                H_(4,2) = 2.0*(BT*BD+B*EV_(3)*CT);
                H_(2,4) = H_(4,2);
                H_(2,2) = 2.0*(BT*BT-B*EV_(3)*EV_(4)*ST);
                varargout{3} = H_;
            end
        end

    case 'eDISTX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ST = sin(EV_(2));
        CT = cos(EV_(2));
        B = EV_(1)+EV_(3)*EV_(4)*CT-pbm.elpar{iel_}(1);
        BA = EV_(4)*CT;
        BD = EV_(3)*CT;
        BT = -EV_(3)*EV_(4)*ST;
        varargout{1} = B*B;
        if(nargout>1)
            g_(1,1) = 2.0*B;
            g_(3,1) = 2.0*B*BA;
            g_(4,1) = 2.0*B*BD;
            g_(2,1) = 2.0*B*BT;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = 2.0;
                H_(1,3) = 2.0*BA;
                H_(3,1) = H_(1,3);
                H_(1,4) = 2.0*BD;
                H_(4,1) = H_(1,4);
                H_(1,2) = 2.0*BT;
                H_(2,1) = H_(1,2);
                H_(3,3) = 2.0*BA*BA;
                H_(3,4) = 2.0*(BD*BA+B*CT);
                H_(4,3) = H_(3,4);
                H_(3,2) = 2.0*(BT*BA-B*EV_(4)*ST);
                H_(2,3) = H_(3,2);
                H_(4,4) = 2.0*BD*BD;
                H_(4,2) = 2.0*(BT*BD-B*EV_(3)*ST);
                H_(2,4) = H_(4,2);
                H_(2,2) = 2.0*(BT*BT-B*EV_(3)*EV_(4)*CT);
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

