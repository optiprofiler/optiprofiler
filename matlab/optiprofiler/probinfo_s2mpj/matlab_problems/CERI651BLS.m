function varargout = CERI651BLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CERI651BLS
%    *********
% 
%    ISIS Data fitting problem CERI651B given as an inconsistent set of
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
%    subset X in [26047.3026604, 26393.719109]
% 
%    SIF input: Nick Gould and Tyrone Rees, Mar 2016
%    Least-squares version of CERI651B.SIF, Nick Gould, Jan 2020.
% 
%    classification = 'SUR2-MN-7-0'
% 
%    Potential and actual number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CERI651BLS';

switch(action)

    case 'setup'

    pb.name      = 'CERI651BLS';
    pb.sifpbname = 'CERI651BLS';
    pbm.name     = 'CERI651BLS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('MPOT') = 10186;
        v_('M') = 66;
        v_('MLOWER') = 7343;
        v_('MUPPER') = 7408;
        v_('N') = 7;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('X7343') = 26052.42188;
        v_('X7344') = 26057.64063;
        v_('X7345') = 26062.85938;
        v_('X7346') = 26068.07813;
        v_('X7347') = 26073.29688;
        v_('X7348') = 26078.51563;
        v_('X7349') = 26083.73438;
        v_('X7350') = 26088.95313;
        v_('X7351') = 26094.17188;
        v_('X7352') = 26099.39063;
        v_('X7353') = 26104.60938;
        v_('X7354') = 26109.82813;
        v_('X7355') = 26115.04688;
        v_('X7356') = 26120.26563;
        v_('X7357') = 26125.48438;
        v_('X7358') = 26130.70313;
        v_('X7359') = 26135.92188;
        v_('X7360') = 26141.14063;
        v_('X7361') = 26146.35938;
        v_('X7362') = 26151.57813;
        v_('X7363') = 26156.79688;
        v_('X7364') = 26162.01563;
        v_('X7365') = 26167.23438;
        v_('X7366') = 26172.45313;
        v_('X7367') = 26177.68750;
        v_('X7368') = 26182.93750;
        v_('X7369') = 26188.18750;
        v_('X7370') = 26193.43750;
        v_('X7371') = 26198.68750;
        v_('X7372') = 26203.93750;
        v_('X7373') = 26209.18750;
        v_('X7374') = 26214.43750;
        v_('X7375') = 26219.68750;
        v_('X7376') = 26224.93750;
        v_('X7377') = 26230.18750;
        v_('X7378') = 26235.43750;
        v_('X7379') = 26240.68750;
        v_('X7380') = 26245.93750;
        v_('X7381') = 26251.18750;
        v_('X7382') = 26256.43750;
        v_('X7383') = 26261.68750;
        v_('X7384') = 26266.93750;
        v_('X7385') = 26272.18750;
        v_('X7386') = 26277.43750;
        v_('X7387') = 26282.68750;
        v_('X7388') = 26287.93750;
        v_('X7389') = 26293.18750;
        v_('X7390') = 26298.43750;
        v_('X7391') = 26303.68750;
        v_('X7392') = 26308.93750;
        v_('X7393') = 26314.18750;
        v_('X7394') = 26319.43750;
        v_('X7395') = 26324.68750;
        v_('X7396') = 26329.93750;
        v_('X7397') = 26335.20313;
        v_('X7398') = 26340.48438;
        v_('X7399') = 26345.76563;
        v_('X7400') = 26351.04688;
        v_('X7401') = 26356.32813;
        v_('X7402') = 26361.60938;
        v_('X7403') = 26366.89063;
        v_('X7404') = 26372.17188;
        v_('X7405') = 26377.45313;
        v_('X7406') = 26382.73438;
        v_('X7407') = 26388.01563;
        v_('X7408') = 26393.29688;
        v_('Y7343') = 1.96083316;
        v_('Y7344') = 0.98041658;
        v_('Y7345') = 0.00000000;
        v_('Y7346') = 1.96083316;
        v_('Y7347') = 3.92166632;
        v_('Y7348') = 0.00000000;
        v_('Y7349') = 3.92166632;
        v_('Y7350') = 0.98041658;
        v_('Y7351') = 1.96083316;
        v_('Y7352') = 1.96083316;
        v_('Y7353') = 3.92166632;
        v_('Y7354') = 1.96083316;
        v_('Y7355') = 0.98041658;
        v_('Y7356') = 4.90208290;
        v_('Y7357') = 8.82374922;
        v_('Y7358') = 5.88249948;
        v_('Y7359') = 14.70624870;
        v_('Y7360') = 12.74541554;
        v_('Y7361') = 27.45166424;
        v_('Y7362') = 27.45166424;
        v_('Y7363') = 32.35374715;
        v_('Y7364') = 52.94249533;
        v_('Y7365') = 43.13832953;
        v_('Y7366') = 47.05999585;
        v_('Y7367') = 50.00124559;
        v_('Y7368') = 61.76624455;
        v_('Y7369') = 52.94249533;
        v_('Y7370') = 34.31458031;
        v_('Y7371') = 42.15791295;
        v_('Y7372') = 32.35374715;
        v_('Y7373') = 40.19707979;
        v_('Y7374') = 33.33416373;
        v_('Y7375') = 23.52999792;
        v_('Y7376') = 16.66708186;
        v_('Y7377') = 12.74541554;
        v_('Y7378') = 13.72583212;
        v_('Y7379') = 13.72583212;
        v_('Y7380') = 12.74541554;
        v_('Y7381') = 6.86291606;
        v_('Y7382') = 14.70624870;
        v_('Y7383') = 7.84333264;
        v_('Y7384') = 8.82374922;
        v_('Y7385') = 8.82374922;
        v_('Y7386') = 6.86291606;
        v_('Y7387') = 9.80416580;
        v_('Y7388') = 4.90208290;
        v_('Y7389') = 4.90208290;
        v_('Y7390') = 3.92166632;
        v_('Y7391') = 1.96083316;
        v_('Y7392') = 2.94124974;
        v_('Y7393') = 2.94124974;
        v_('Y7394') = 0.00000000;
        v_('Y7395') = 0.98041658;
        v_('Y7396') = 0.98041658;
        v_('Y7397') = 3.92166632;
        v_('Y7398') = 0.98041658;
        v_('Y7399') = 0.98041658;
        v_('Y7400') = 1.96083316;
        v_('Y7401') = 1.96083316;
        v_('Y7402') = 0.00000000;
        v_('Y7403') = 3.92166632;
        v_('Y7404') = 0.98041658;
        v_('Y7405') = 0.98041658;
        v_('Y7406') = 0.98041658;
        v_('Y7407') = 2.94124974;
        v_('Y7408') = 0.00000000;
        v_('E7343') = 1.41421356;
        v_('E7344') = 1.00000000;
        v_('E7345') = 1.00000000;
        v_('E7346') = 1.41421356;
        v_('E7347') = 2.00000000;
        v_('E7348') = 1.00000000;
        v_('E7349') = 2.00000000;
        v_('E7350') = 1.00000000;
        v_('E7351') = 1.41421356;
        v_('E7352') = 1.41421356;
        v_('E7353') = 2.00000000;
        v_('E7354') = 1.41421356;
        v_('E7355') = 1.00000000;
        v_('E7356') = 2.23606798;
        v_('E7357') = 3.00000000;
        v_('E7358') = 2.44948974;
        v_('E7359') = 3.87298335;
        v_('E7360') = 3.60555128;
        v_('E7361') = 5.29150262;
        v_('E7362') = 5.29150262;
        v_('E7363') = 5.74456265;
        v_('E7364') = 7.34846923;
        v_('E7365') = 6.63324958;
        v_('E7366') = 6.92820323;
        v_('E7367') = 7.14142843;
        v_('E7368') = 7.93725393;
        v_('E7369') = 7.34846923;
        v_('E7370') = 5.91607978;
        v_('E7371') = 6.55743852;
        v_('E7372') = 5.74456265;
        v_('E7373') = 6.40312424;
        v_('E7374') = 5.83095189;
        v_('E7375') = 4.89897949;
        v_('E7376') = 4.12310563;
        v_('E7377') = 3.60555128;
        v_('E7378') = 3.74165739;
        v_('E7379') = 3.74165739;
        v_('E7380') = 3.60555128;
        v_('E7381') = 2.64575131;
        v_('E7382') = 3.87298335;
        v_('E7383') = 2.82842712;
        v_('E7384') = 3.00000000;
        v_('E7385') = 3.00000000;
        v_('E7386') = 2.64575131;
        v_('E7387') = 3.16227766;
        v_('E7388') = 2.23606798;
        v_('E7389') = 2.23606798;
        v_('E7390') = 2.00000000;
        v_('E7391') = 1.41421356;
        v_('E7392') = 1.73205081;
        v_('E7393') = 1.73205081;
        v_('E7394') = 1.00000000;
        v_('E7395') = 1.00000000;
        v_('E7396') = 1.00000000;
        v_('E7397') = 2.00000000;
        v_('E7398') = 1.00000000;
        v_('E7399') = 1.00000000;
        v_('E7400') = 1.41421356;
        v_('E7401') = 1.41421356;
        v_('E7402') = 1.00000000;
        v_('E7403') = 2.00000000;
        v_('E7404') = 1.00000000;
        v_('E7405') = 1.00000000;
        v_('E7406') = 1.00000000;
        v_('E7407') = 1.73205081;
        v_('E7408') = 1.00000000;
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
        pb.x0(ix_('I'),1) = 3527.31;
        pb.x0(ix_('S'),1) = 29.4219;
        pb.x0(ix_('X0'),1) = 26185.9;
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

