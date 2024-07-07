function varargout = CERI651CLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CERI651CLS
%    *********
% 
%    ISIS Data fitting problem CERI651C given as an inconsistent set of
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
%    subset X in [23919.5789114, 24189.3183142]
% 
%    SIF input: Nick Gould and Tyrone Rees, Mar 2016
%    Least-squares version of CERI651C.SIF, Nick Gould, Jan 2020.
% 
%    classification = 'SUR2-MN-7-0'
% 
%    Potential and actual number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CERI651CLS';

switch(action)

    case 'setup'

    pb.name      = 'CERI651CLS';
    pb.sifpbname = 'CERI651CLS';
    pbm.name     = 'CERI651CLS';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('MPOT') = 10186;
        v_('M') = 56;
        v_('MLOWER') = 6916;
        v_('MUPPER') = 6971;
        v_('N') = 7;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('X6916') = 23920.10938;
        v_('X6917') = 23924.89063;
        v_('X6918') = 23929.67188;
        v_('X6919') = 23934.45313;
        v_('X6920') = 23939.23438;
        v_('X6921') = 23944.01563;
        v_('X6922') = 23948.79688;
        v_('X6923') = 23953.57813;
        v_('X6924') = 23958.35938;
        v_('X6925') = 23963.14063;
        v_('X6926') = 23967.92188;
        v_('X6927') = 23972.70313;
        v_('X6928') = 23977.48438;
        v_('X6929') = 23982.26563;
        v_('X6930') = 23987.06250;
        v_('X6931') = 23991.87500;
        v_('X6932') = 23996.68750;
        v_('X6933') = 24001.50000;
        v_('X6934') = 24006.31250;
        v_('X6935') = 24011.12500;
        v_('X6936') = 24015.93750;
        v_('X6937') = 24020.75000;
        v_('X6938') = 24025.56250;
        v_('X6939') = 24030.37500;
        v_('X6940') = 24035.18750;
        v_('X6941') = 24040.00000;
        v_('X6942') = 24044.81250;
        v_('X6943') = 24049.62500;
        v_('X6944') = 24054.43750;
        v_('X6945') = 24059.25000;
        v_('X6946') = 24064.06250;
        v_('X6947') = 24068.87500;
        v_('X6948') = 24073.68750;
        v_('X6949') = 24078.50000;
        v_('X6950') = 24083.31250;
        v_('X6951') = 24088.12500;
        v_('X6952') = 24092.93750;
        v_('X6953') = 24097.75000;
        v_('X6954') = 24102.56250;
        v_('X6955') = 24107.37500;
        v_('X6956') = 24112.18750;
        v_('X6957') = 24117.00000;
        v_('X6958') = 24121.81250;
        v_('X6959') = 24126.62500;
        v_('X6960') = 24131.43750;
        v_('X6961') = 24136.25000;
        v_('X6962') = 24141.06250;
        v_('X6963') = 24145.89063;
        v_('X6964') = 24150.73438;
        v_('X6965') = 24155.57813;
        v_('X6966') = 24160.42188;
        v_('X6967') = 24165.26563;
        v_('X6968') = 24170.10938;
        v_('X6969') = 24174.95313;
        v_('X6970') = 24179.79688;
        v_('X6971') = 24184.64063;
        v_('Y6916') = 0.00000000;
        v_('Y6917') = 0.98041658;
        v_('Y6918') = 1.96083316;
        v_('Y6919') = 0.00000000;
        v_('Y6920') = 0.98041658;
        v_('Y6921') = 0.00000000;
        v_('Y6922') = 0.00000000;
        v_('Y6923') = 3.92166632;
        v_('Y6924') = 0.98041658;
        v_('Y6925') = 0.00000000;
        v_('Y6926') = 0.98041658;
        v_('Y6927') = 2.94124974;
        v_('Y6928') = 1.96083316;
        v_('Y6929') = 0.98041658;
        v_('Y6930') = 2.94124974;
        v_('Y6931') = 8.82374922;
        v_('Y6932') = 5.88249948;
        v_('Y6933') = 6.86291606;
        v_('Y6934') = 8.82374922;
        v_('Y6935') = 11.76499896;
        v_('Y6936') = 12.74541554;
        v_('Y6937') = 6.86291606;
        v_('Y6938') = 8.82374922;
        v_('Y6939') = 12.74541554;
        v_('Y6940') = 13.72583212;
        v_('Y6941') = 8.82374922;
        v_('Y6942') = 12.74541554;
        v_('Y6943') = 19.60833160;
        v_('Y6944') = 4.90208290;
        v_('Y6945') = 2.94124974;
        v_('Y6946') = 1.96083316;
        v_('Y6947') = 3.92166632;
        v_('Y6948') = 3.92166632;
        v_('Y6949') = 5.88249948;
        v_('Y6950') = 2.94124974;
        v_('Y6951') = 4.90208290;
        v_('Y6952') = 6.86291606;
        v_('Y6953') = 2.94124974;
        v_('Y6954') = 1.96083316;
        v_('Y6955') = 0.00000000;
        v_('Y6956') = 1.96083316;
        v_('Y6957') = 2.94124974;
        v_('Y6958') = 1.96083316;
        v_('Y6959') = 1.96083316;
        v_('Y6960') = 1.96083316;
        v_('Y6961') = 3.92166632;
        v_('Y6962') = 0.00000000;
        v_('Y6963') = 0.00000000;
        v_('Y6964') = 3.92166632;
        v_('Y6965') = 2.94124974;
        v_('Y6966') = 1.96083316;
        v_('Y6967') = 0.00000000;
        v_('Y6968') = 1.96083316;
        v_('Y6969') = 0.00000000;
        v_('Y6970') = 0.98041658;
        v_('Y6971') = 0.98041658;
        v_('E6916') = 1.00000000;
        v_('E6917') = 1.00000000;
        v_('E6918') = 1.41421356;
        v_('E6919') = 1.00000000;
        v_('E6920') = 1.00000000;
        v_('E6921') = 1.00000000;
        v_('E6922') = 1.00000000;
        v_('E6923') = 2.00000000;
        v_('E6924') = 1.00000000;
        v_('E6925') = 1.00000000;
        v_('E6926') = 1.00000000;
        v_('E6927') = 1.73205081;
        v_('E6928') = 1.41421356;
        v_('E6929') = 1.00000000;
        v_('E6930') = 1.73205081;
        v_('E6931') = 3.00000000;
        v_('E6932') = 2.44948974;
        v_('E6933') = 2.64575131;
        v_('E6934') = 3.00000000;
        v_('E6935') = 3.46410162;
        v_('E6936') = 3.60555128;
        v_('E6937') = 2.64575131;
        v_('E6938') = 3.00000000;
        v_('E6939') = 3.60555128;
        v_('E6940') = 3.74165739;
        v_('E6941') = 3.00000000;
        v_('E6942') = 3.60555128;
        v_('E6943') = 4.47213595;
        v_('E6944') = 2.23606798;
        v_('E6945') = 1.73205081;
        v_('E6946') = 1.41421356;
        v_('E6947') = 2.00000000;
        v_('E6948') = 2.00000000;
        v_('E6949') = 2.44948974;
        v_('E6950') = 1.73205081;
        v_('E6951') = 2.23606798;
        v_('E6952') = 2.64575131;
        v_('E6953') = 1.73205081;
        v_('E6954') = 1.41421356;
        v_('E6955') = 1.00000000;
        v_('E6956') = 1.41421356;
        v_('E6957') = 1.73205081;
        v_('E6958') = 1.41421356;
        v_('E6959') = 1.41421356;
        v_('E6960') = 1.41421356;
        v_('E6961') = 2.00000000;
        v_('E6962') = 1.00000000;
        v_('E6963') = 1.00000000;
        v_('E6964') = 2.00000000;
        v_('E6965') = 1.73205081;
        v_('E6966') = 1.41421356;
        v_('E6967') = 1.00000000;
        v_('E6968') = 1.41421356;
        v_('E6969') = 1.00000000;
        v_('E6970') = 1.00000000;
        v_('E6971') = 1.00000000;
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
        pb.x0(ix_('I'),1) = 597.076;
        pb.x0(ix_('S'),1) = 22.9096;
        pb.x0(ix_('X0'),1) = 24027.5;
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

