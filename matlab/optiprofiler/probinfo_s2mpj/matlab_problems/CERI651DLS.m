function varargout = CERI651DLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CERI651DLS
%    *********
% 
%    ISIS Data fitting problem CERI651D given as an inconsistent set of
%    nonlinear equations.
% 
%    Fit: y = c + l * x + I*A*B/2(A+B) *
%               [ exp( A*[A*S^2+2(x-X0)]/2) * erfc( A*S^2+(x-X0)/S*sqrt(2) ) +
%                 exp( B*[B*S^2+2(x-X0)]/2) * erfc( B*S^2+(x-X0)/S*sqrt(2) ) ]
% 
%    Source: fit to a sum of a linear background and a back-to-back exponential
%    using data enginx_ceria193749_spectrum_number_651_vana_corrected-0
%    from Mantid (http://www.mantidproject.org)
% 
%    subset X in [12986.356148, 13161.356148]
% 
%    SIF input: Nick Gould and Tyrone Rees, Mar 2016
%    Least-squares version of CERI651D.SIF, Nick Gould, Jan 2020.
% 
%    classification = 'SUR2-MN-7-0'
% 
%    Potential and actual number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CERI651DLS';

switch(action)

    case 'setup'

    pb.name      = 'CERI651DLS';
    pb.sifpbname = 'CERI651DLS';
    pbm.name     = 'CERI651DLS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('MPOT') = 10186;
        v_('M') = 67;
        v_('MLOWER') = 3862;
        v_('MUPPER') = 3928;
        v_('N') = 7;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('X3862') = 12987.48438;
        v_('X3863') = 12990.07813;
        v_('X3864') = 12992.67188;
        v_('X3865') = 12995.26563;
        v_('X3866') = 12997.85938;
        v_('X3867') = 13000.45313;
        v_('X3868') = 13003.04688;
        v_('X3869') = 13005.64063;
        v_('X3870') = 13008.23438;
        v_('X3871') = 13010.82813;
        v_('X3872') = 13013.42188;
        v_('X3873') = 13016.01563;
        v_('X3874') = 13018.60938;
        v_('X3875') = 13021.20313;
        v_('X3876') = 13023.79688;
        v_('X3877') = 13026.39063;
        v_('X3878') = 13028.98438;
        v_('X3879') = 13031.57813;
        v_('X3880') = 13034.17188;
        v_('X3881') = 13036.76563;
        v_('X3882') = 13039.35938;
        v_('X3883') = 13041.95313;
        v_('X3884') = 13044.54688;
        v_('X3885') = 13047.14063;
        v_('X3886') = 13049.75000;
        v_('X3887') = 13052.37500;
        v_('X3888') = 13055.00000;
        v_('X3889') = 13057.62500;
        v_('X3890') = 13060.25000;
        v_('X3891') = 13062.87500;
        v_('X3892') = 13065.50000;
        v_('X3893') = 13068.12500;
        v_('X3894') = 13070.75000;
        v_('X3895') = 13073.37500;
        v_('X3896') = 13076.00000;
        v_('X3897') = 13078.62500;
        v_('X3898') = 13081.25000;
        v_('X3899') = 13083.87500;
        v_('X3900') = 13086.50000;
        v_('X3901') = 13089.12500;
        v_('X3902') = 13091.75000;
        v_('X3903') = 13094.37500;
        v_('X3904') = 13097.00000;
        v_('X3905') = 13099.62500;
        v_('X3906') = 13102.25000;
        v_('X3907') = 13104.87500;
        v_('X3908') = 13107.50000;
        v_('X3909') = 13110.12500;
        v_('X3910') = 13112.75000;
        v_('X3911') = 13115.37500;
        v_('X3912') = 13118.00000;
        v_('X3913') = 13120.62500;
        v_('X3914') = 13123.25000;
        v_('X3915') = 13125.87500;
        v_('X3916') = 13128.50000;
        v_('X3917') = 13131.12500;
        v_('X3918') = 13133.75000;
        v_('X3919') = 13136.37500;
        v_('X3920') = 13139.00000;
        v_('X3921') = 13141.62500;
        v_('X3922') = 13144.25000;
        v_('X3923') = 13146.87500;
        v_('X3924') = 13149.50000;
        v_('X3925') = 13152.12500;
        v_('X3926') = 13154.75000;
        v_('X3927') = 13157.37500;
        v_('X3928') = 13160.00000;
        v_('Y3862') = 0.00000000;
        v_('Y3863') = 0.00000000;
        v_('Y3864') = 0.00000000;
        v_('Y3865') = 0.00000000;
        v_('Y3866') = 0.00000000;
        v_('Y3867') = 0.00000000;
        v_('Y3868') = 0.00000000;
        v_('Y3869') = 1.96083316;
        v_('Y3870') = 0.00000000;
        v_('Y3871') = 0.00000000;
        v_('Y3872') = 0.00000000;
        v_('Y3873') = 0.00000000;
        v_('Y3874') = 0.00000000;
        v_('Y3875') = 0.00000000;
        v_('Y3876') = 0.00000000;
        v_('Y3877') = 0.00000000;
        v_('Y3878') = 0.00000000;
        v_('Y3879') = 0.00000000;
        v_('Y3880') = 0.00000000;
        v_('Y3881') = 0.00000000;
        v_('Y3882') = 0.00000000;
        v_('Y3883') = 0.00000000;
        v_('Y3884') = 0.00000000;
        v_('Y3885') = 0.00000000;
        v_('Y3886') = 0.00000000;
        v_('Y3887') = 0.98041658;
        v_('Y3888') = 0.00000000;
        v_('Y3889') = 0.00000000;
        v_('Y3890') = 0.00000000;
        v_('Y3891') = 0.00000000;
        v_('Y3892') = 0.00000000;
        v_('Y3893') = 0.00000000;
        v_('Y3894') = 0.00000000;
        v_('Y3895') = 4.90208290;
        v_('Y3896') = 0.98041658;
        v_('Y3897') = 0.98041658;
        v_('Y3898') = 0.98041658;
        v_('Y3899') = 3.92166632;
        v_('Y3900') = 1.96083316;
        v_('Y3901') = 1.96083316;
        v_('Y3902') = 0.98041658;
        v_('Y3903') = 1.96083316;
        v_('Y3904') = 1.96083316;
        v_('Y3905') = 1.96083316;
        v_('Y3906') = 0.98041658;
        v_('Y3907') = 0.98041658;
        v_('Y3908') = 0.00000000;
        v_('Y3909') = 0.00000000;
        v_('Y3910') = 0.00000000;
        v_('Y3911') = 0.00000000;
        v_('Y3912') = 0.98041658;
        v_('Y3913') = 0.98041658;
        v_('Y3914') = 0.00000000;
        v_('Y3915') = 0.00000000;
        v_('Y3916') = 0.00000000;
        v_('Y3917') = 0.00000000;
        v_('Y3918') = 0.00000000;
        v_('Y3919') = 0.00000000;
        v_('Y3920') = 1.96083316;
        v_('Y3921') = 0.98041658;
        v_('Y3922') = 0.00000000;
        v_('Y3923') = 0.00000000;
        v_('Y3924') = 0.00000000;
        v_('Y3925') = 0.00000000;
        v_('Y3926') = 0.00000000;
        v_('Y3927') = 0.00000000;
        v_('Y3928') = 0.00000000;
        v_('E3862') = 1.00000000;
        v_('E3863') = 1.00000000;
        v_('E3864') = 1.00000000;
        v_('E3865') = 1.00000000;
        v_('E3866') = 1.00000000;
        v_('E3867') = 1.00000000;
        v_('E3868') = 1.00000000;
        v_('E3869') = 1.41421356;
        v_('E3870') = 1.00000000;
        v_('E3871') = 1.00000000;
        v_('E3872') = 1.00000000;
        v_('E3873') = 1.00000000;
        v_('E3874') = 1.00000000;
        v_('E3875') = 1.00000000;
        v_('E3876') = 1.00000000;
        v_('E3877') = 1.00000000;
        v_('E3878') = 1.00000000;
        v_('E3879') = 1.00000000;
        v_('E3880') = 1.00000000;
        v_('E3881') = 1.00000000;
        v_('E3882') = 1.00000000;
        v_('E3883') = 1.00000000;
        v_('E3884') = 1.00000000;
        v_('E3885') = 1.00000000;
        v_('E3886') = 1.00000000;
        v_('E3887') = 1.00000000;
        v_('E3888') = 1.00000000;
        v_('E3889') = 1.00000000;
        v_('E3890') = 1.00000000;
        v_('E3891') = 1.00000000;
        v_('E3892') = 1.00000000;
        v_('E3893') = 1.00000000;
        v_('E3894') = 1.00000000;
        v_('E3895') = 2.23606798;
        v_('E3896') = 1.00000000;
        v_('E3897') = 1.00000000;
        v_('E3898') = 1.00000000;
        v_('E3899') = 2.00000000;
        v_('E3900') = 1.41421356;
        v_('E3901') = 1.41421356;
        v_('E3902') = 1.00000000;
        v_('E3903') = 1.41421356;
        v_('E3904') = 1.41421356;
        v_('E3905') = 1.41421356;
        v_('E3906') = 1.00000000;
        v_('E3907') = 1.00000000;
        v_('E3908') = 1.00000000;
        v_('E3909') = 1.00000000;
        v_('E3910') = 1.00000000;
        v_('E3911') = 1.00000000;
        v_('E3912') = 1.00000000;
        v_('E3913') = 1.00000000;
        v_('E3914') = 1.00000000;
        v_('E3915') = 1.00000000;
        v_('E3916') = 1.00000000;
        v_('E3917') = 1.00000000;
        v_('E3918') = 1.00000000;
        v_('E3919') = 1.00000000;
        v_('E3920') = 1.41421356;
        v_('E3921') = 1.00000000;
        v_('E3922') = 1.00000000;
        v_('E3923') = 1.00000000;
        v_('E3924') = 1.00000000;
        v_('E3925') = 1.00000000;
        v_('E3926') = 1.00000000;
        v_('E3927') = 1.00000000;
        v_('E3928') = 1.00000000;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2xlib('ii','C',ix_);
        pb.xnames{iv} = 'C';
        [iv,ix_] = s2xlib('ii','L',ix_);
        pb.xnames{iv} = 'L';
        [iv,ix_] = s2xlib('ii','A',ix_);
        pb.xnames{iv} = 'A';
        [iv,ix_] = s2xlib('ii','B',ix_);
        pb.xnames{iv} = 'B';
        [iv,ix_] = s2xlib('ii','I',ix_);
        pb.xnames{iv} = 'I';
        [iv,ix_] = s2xlib('ii','S',ix_);
        pb.xnames{iv} = 'S';
        [iv,ix_] = s2xlib('ii','X0',ix_);
        pb.xnames{iv} = 'X0';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('MLOWER'):v_('MUPPER')
            v_('E') = v_(['E',int2str(I)]);
            v_('EINV') = v_('ONE')/v_('E');
            v_('XOVERE') = v_('EINV')*v_(['X',int2str(I)]);
            [ig,ig_] = s2xlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('C');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('EINV')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('EINV');
            end
            iv = ix_('L');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XOVERE')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XOVERE');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('MLOWER'):v_('MUPPER')
            v_('E') = v_(['E',int2str(I)]);
            v_('EINV') = v_('ONE')/v_('E');
            v_('YOVERE') = v_('EINV')*v_(['Y',int2str(I)]);
            pbm.gconst(ig_(['F',int2str(I)])) = v_('YOVERE');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('C'),1) = 0.0;
        pb.x0(ix_('L'),1) = 0.0;
        pb.x0(ix_('A'),1) = 1.0;
        pb.x0(ix_('B'),1) = 0.05;
        pb.x0(ix_('I'),1) = 15.1595;
        pb.x0(ix_('S'),1) = 8.0;
        pb.x0(ix_('X0'),1) = 13072.9;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'B2BEXP',iet_);
        elftv{it}{1} = 'A';
        elftv{it}{2} = 'B';
        elftv{it}{3} = 'I';
        elftv{it}{4} = 'S';
        elftv{it}{5} = 'Y';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('MLOWER'):v_('MUPPER')
            ename = ['B',int2str(I)];
            [ie,ie_] = s2xlib('ii',ename,ie_);
            pbm.elftype{ie} = 'B2BEXP';
            ielftype(ie) = iet_('B2BEXP');
            vname = 'A';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('A',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = I;
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('I',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'S';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('S',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X0';
            [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','L2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'L2';
        end
        for I=v_('MLOWER'):v_('MUPPER')
            v_('E') = v_(['E',int2str(I)]);
            v_('EINV') = v_('ONE')/v_('E');
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = v_('EINV');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        pb.pbclass = 'SUR2-MN-7-0';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = sqrt(1.0e0/atan(1.0e0));    % this is TORPI
        pbm.efpar(2) = sqrt(0.5e0);    % this is ROOTP5
        varargout{1} = pbm;

    case 'B2BEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        APB = EV_(1)+EV_(2);
        APB2 = APB*APB;
        A2 = EV_(1)*EV_(1);
        B2 = EV_(2)*EV_(2);
        AB = EV_(1)*EV_(2);
        S2 = EV_(4)*EV_(4);
        S3 = EV_(4)*S2;
        PI = 0.5e0*AB/APB;
        P = EV_(3)*PI;
        PAI = 0.5e0*EV_(2)/APB-0.5e0*AB/APB2;
        PBI = 0.5e0*EV_(1)/APB-0.5e0*AB/APB2;
        PA = PAI*EV_(3);
        loc_PB = PBI*EV_(3);
        PAB = EV_(3)*AB/APB^3;
        PAA = -EV_(3)*EV_(2)/APB2+PAB;
        PBB = -EV_(3)*EV_(1)/APB2+PAB;
        XMY = pbm.elpar{iel_}(1)-EV_(5);
        Z = XMY/EV_(4);
        ZY = -1.0e0/EV_(4);
        ZS = -XMY/EV_(4)^2;
        ZSY = 1.0e0/EV_(4)^2;
        ZSS = 2.0e0*XMY/EV_(4)^3;
        R = exp(-0.5e0*Z^2);
        DR = -Z*R;
        D2R = -R-Z*DR;
        RS = DR*ZS;
        RY = DR*ZY;
        RSS = D2R*ZS*ZS+DR*ZSS;
        RSY = D2R*ZS*ZY+DR*ZSY;
        RYY = D2R*ZY*ZY;
        AC = pbm.efpar(2)*(EV_(1)*EV_(4)+XMY/EV_(4));
        ACA = pbm.efpar(2)*EV_(4);
        ACS = pbm.efpar(2)*(EV_(1)-XMY/S2);
        ACY = -pbm.efpar(2)/EV_(4);
        ACAS = pbm.efpar(2);
        ACSS = 2.0e0*pbm.efpar(2)*XMY/S3;
        ACSY = pbm.efpar(2)/S2;
        BC = pbm.efpar(2)*(EV_(2)*EV_(4)+XMY/EV_(4));
        BCB = ACA;
        BCS = pbm.efpar(2)*(EV_(2)-XMY/S2);
        BCY = ACY;
        BCBS = pbm.efpar(2);
        BCSS = ACSS;
        BCSY = ACSY;
        QA = ERFC_EV_(4)CEV_(1)LED(AC);
        DQA = 2.0e0*AC*QA-pbm.efpar(1);
        D2QA = 2.0e0*(QA+AC*DQA);
        QAA = DQA*ACA;
        QAS = DQA*ACS;
        QAY = DQA*ACY;
        QAAA = D2QA*ACA*ACA;
        QAAS = D2QA*ACA*ACS+DQA*ACAS;
        QAAY = D2QA*ACA*ACY;
        QASS = D2QA*ACS*ACS+DQA*ACSS;
        QASY = D2QA*ACS*ACY+DQA*ACSY;
        QAYY = D2QA*ACY*ACY;
        QB = ERFC_EV_(4)CEV_(1)LED(BC);
        DQB = 2.0e0*BC*QB-pbm.efpar(1);
        D2QB = 2.0e0*(QB+BC*DQB);
        QBB = DQB*BCB;
        QBS = DQB*BCS;
        QBY = DQB*BCY;
        QBBB = D2QB*BCB*BCB;
        QBBS = D2QB*BCB*BCS+DQB*BCBS;
        QBBY = D2QB*BCB*BCY;
        QBSS = D2QB*BCS*BCS+DQB*BCSS;
        QBSY = D2QB*BCS*BCY+DQB*BCSY;
        QBYY = D2QB*BCY*BCY;
        T = QA+QB;
        TA = QAA;
        TB = QBB;
        TS = QAS+QBS;
        TY = QAY+QBY;
        TAA = QAAA;
        TAS = QAAS;
        TAY = QAAY;
        TBB = QBBB;
        TBS = QBBS;
        TBY = QBBY;
        TSS = QASS+QBSS;
        TSY = QASY+QBSY;
        TYY = QAYY+QBYY;
        varargout{1} = P*T*R;
        if(nargout>1)
            g_(1,1) = (P*TA+PA*T)*R;
            g_(2,1) = (P*TB+loc_PB*T)*R;
            g_(3,1) = PI*T*R;
            g_(4,1) = P*(T*RS+TS*R);
            g_(5,1) = P*(T*RY+TY*R);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = R*(P*TAA+PAA*T+2.0e0*PA*TA);
                H_(1,2) = R*(PA*TB+loc_PB*TA+PAB*T);
                H_(2,1) = H_(1,2);
                H_(1,3) = (PI*TA+PAI*T)*R;
                H_(3,1) = H_(1,3);
                H_(1,4) = (P*TA+PA*T)*RS+(P*TAS+PA*TS)*R;
                H_(4,1) = H_(1,4);
                H_(1,5) = (P*TA+PA*T)*RY+(P*TAY+PA*TY)*R;
                H_(5,1) = H_(1,5);
                H_(2,2) = R*(P*TBB+PBB*T+2.0e0*loc_PB*TB);
                H_(2,3) = (PI*TB+PBI*T)*R;
                H_(3,2) = H_(2,3);
                H_(2,4) = (P*TB+loc_PB*T)*RS+(P*TBS+loc_PB*TS)*R;
                H_(4,2) = H_(2,4);
                H_(2,5) = (P*TB+loc_PB*T)*RY+(P*TBY+loc_PB*TY)*R;
                H_(5,2) = H_(2,5);
                H_(3,4) = PI*(T*RS+TS*R);
                H_(4,3) = H_(3,4);
                H_(3,5) = PI*(T*RY+TY*R);
                H_(5,3) = H_(3,5);
                H_(4,4) = P*(T*RSS+TSS*R+2.0e0*TS*RS);
                H_(4,5) = P*(T*RSY+TSY*R+TS*RY+TY*RS);
                H_(5,4) = H_(4,5);
                H_(5,5) = P*(T*RYY+TYY*R+2.0e0*TY*RY);
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
                H_ = 2.0e0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [2,0];
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

