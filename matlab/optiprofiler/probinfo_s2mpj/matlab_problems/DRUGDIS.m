function varargout = DRUGDIS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DRUGDIS
%    *********
% 
%    A control problem based on the kinetic model of Aarons and Rowland for
%    DRUG DISplacemnt, which simulates the interaction of the two drugs 
%    (warfarin and phenylnutazone) in a patient bloodstream.  
%    The state variable are the concentrations of unbound warfarin (w) and 
%    phenylbutazone (p).  The problem is to control the rate of injection (u) 
%    of the pain-killing phenylbutazone so that both drugs reach a specified 
%    steady-state in minimum time and the concentration of warfarin does not 
%    rise above a given toxicity level.
% 
%    The problem is discretized using the trapezoidal rule.  It is non-convex.
% 
%    The problem can be made harder by diminishing the value of the lower bound
%    on the final time TF (while maintaining it strictly positive).
% 
%    Source:
%    H. Maurer and M. Wiegand,
%    "Numerical solution of a drug displacement problem with bounded state
%    variables",
%    Optimal Control Applications and Methods 13, pp. 43-55, 1992.
% 
%    SIF input: Ph. Toint, Nov 1993.
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'LOR2-MN-V-V'
% 
%    Discretization: specify the number of interior points + 1
% 
%       Alternative values for the SIF file parameters:
% IE NI                  10             $-PARAMETER n=  34, m= 20 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DRUGDIS';

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
        if(nargs<1)
            v_('NI') = 10;  %  SIF file default value
        else
            v_('NI') = varargin{1};
        end
% IE NI                  50             $-PARAMETER n= 154, m=100 
% IE NI                  100            $-PARAMETER n= 304, m=200  original value
% IE NI                  200            $-PARAMETER n= 604, m=400 
% IE NI                  500            $-PARAMETER n=1504, m=1000 
% IE NI                  1000           $-PARAMETER n=3004, m=2000 
% IE NI                  2000           $-PARAMETER n=6004, m=4000 
        if(nargs<2)
            v_('TOXIC') = 0.026;  %  SIF file default value
        else
            v_('TOXIC') = varargin{2};
        end
        if(nargs<3)
            v_('WSS') = 0.02;  %  SIF file default value
        else
            v_('WSS') = varargin{3};
        end
        if(nargs<4)
            v_('UMAX') = 8.0;  %  SIF file default value
        else
            v_('UMAX') = varargin{4};
        end
        if(nargs<5)
            v_('PSTART') = 0.0;  %  SIF file default value
        else
            v_('PSTART') = varargin{5};
        end
        if(nargs<6)
            v_('PFINAL') = 2.0;  %  SIF file default value
        else
            v_('PFINAL') = varargin{6};
        end
        v_('AVP') = v_('PSTART')+v_('PFINAL');
        v_('AVP') = 0.5*v_('AVP');
        v_('NI-1') = -1+v_('NI');
        v_('RNI') = v_('NI');
        v_('-1/2NI') = -0.5/v_('RNI');
        v_('0') = 0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','TF',ix_);
        pb.xnames{iv} = 'TF';
        pb.xscale(iv,1) = 200.0;
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I)],ix_);
            pb.xnames{iv} = ['W',int2str(I)];
            pb.xscale(iv,1) = 0.02;
        end
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for I=v_('0'):v_('NI')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','TFINAL',ig_);
        gtype{ig} = '<>';
        iv = ix_('TF');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        pbm.gscale(ig,1) = 100.0;
        for I=v_('0'):v_('NI-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['EW',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EW',int2str(I)];
            iv = ix_(['W',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['W',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            pbm.gscale(ig,1) = 0.02;
            [ig,ig_] = s2mpjlib('ii',['EP',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EP',int2str(I)];
            iv = ix_(['P',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['P',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('TF'),1) = 200.0;
        for I=v_('0'):v_('NI')
            pb.xupper(ix_(['W',int2str(I)])) = v_('TOXIC');
        end
        for I=v_('0'):v_('NI-1')
            pb.xupper(ix_(['U',int2str(I)])) = v_('UMAX');
        end
        pb.xlower(ix_(['W',int2str(round(v_('0')))]),1) = v_('WSS');
        pb.xupper(ix_(['W',int2str(round(v_('0')))]),1) = v_('WSS');
        pb.xlower(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.xupper(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.xlower(ix_(['P',int2str(round(v_('0')))]),1) = v_('PSTART');
        pb.xupper(ix_(['P',int2str(round(v_('0')))]),1) = v_('PSTART');
        pb.xlower(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        pb.xupper(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        v_('DP') = v_('PFINAL')-v_('PSTART');
        v_('DP/NI') = v_('DP')/v_('RNI');
        for I=v_('0'):v_('NI-1')
            v_('RI') = I;
            v_('IDP/NI') = v_('RI')*v_('DP/NI');
            pb.x0(ix_(['P',int2str(I)]),1) = v_('IDP/NI');
            pb.x0(ix_(['W',int2str(I)]),1) = v_('WSS');
            pb.x0(ix_(['U',int2str(I)]),1) = v_('UMAX');
        end
        pb.x0(ix_('TF'),1) = 240.0;
        pb.x0(ix_(['W',int2str(round(v_('NI')))]),1) = v_('WSS');
        pb.x0(ix_(['P',int2str(round(v_('NI')))]),1) = v_('PFINAL');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEW',iet_);
        elftv{it}{1} = 'T';
        elftv{it}{2} = 'W';
        elftv{it}{3} = 'P';
        elftv{it}{4} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'eEP',iet_);
        elftv{it}{1} = 'T';
        elftv{it}{2} = 'W';
        elftv{it}{3} = 'P';
        elftv{it}{4} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('NI')
            ename = ['WA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEW';
            ielftype(ie) = iet_('eEW');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('P',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['PA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEP';
            ielftype(ie) = iet_('eEP');
            vname = 'TF';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['W',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('P',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('0'):v_('NI-1')
            ig = ig_(['EW',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WA',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/2NI');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['WA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/2NI');
            ig = ig_(['EP',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PA',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/2NI');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/2NI');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 200.0;
%    Solution
% LO SOLTN(10)           3.82432
% LO SOLTN(50)           4.19953
% LO SOLTN(100)          4.23934
% LO SOLTN(200)          4.25762
% LO SOLTN(500)
% LO SOLTN(1000)
% LO SOLTN(Maurer)       2.62637
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LOR2-MN-V-V';
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 46.4;
        pbm.efpar(2) = 0.02;
        pbm.efpar(3) = 0.2;
        pbm.efpar(4) = 232.0;
        pbm.efpar(5) = pbm.efpar(1)*pbm.efpar(1);
        varargout{1} = pbm;

    case 'eEW'

        EV_  = varargin{1};
        iel_ = varargin{2};
        D = 1.0+pbm.efpar(3)*(EV_(2)+EV_(3));
        DD = D*D;
        DD1 = 2.0*pbm.efpar(3)*D;
        DD2 = 2.0*pbm.efpar(3)*pbm.efpar(3);
        A = DD+pbm.efpar(4)+pbm.efpar(1)*EV_(2);
        AW = DD1+pbm.efpar(1);
        B = DD+pbm.efpar(4)+pbm.efpar(1)*EV_(3);
        BP = DD1+pbm.efpar(1);
        C = A*B-pbm.efpar(5)*EV_(2)*EV_(3);
        CW = AW*B+A*DD1-pbm.efpar(5)*EV_(3);
        CP = DD1*B+A*BP-pbm.efpar(5)*EV_(2);
        CWW = DD2*B+2.0*AW*DD1+A*DD2;
        CWP = DD2*B+AW*BP+DD1*DD1+A*DD2-pbm.efpar(5);
        CPP = DD2*B+2.0*DD1*BP+A*DD2;
        F = DD/C;
        H = DD1-F*CW;
        I = DD1-F*CP;
        FW = H/C;
        FP = I/C;
        HW = DD2-CW*FW-F*CWW;
        HP = DD2-CW*FW-F*CWP;
        IP = DD2-CP*FP-F*CPP;
        FWW = (HW-FW*CW)/C;
        FWP = (HP-FW*CP)/C;
        FPP = (IP-FP*CP)/C;
        GU = pbm.efpar(1)*EV_(2);
        G = A*(pbm.efpar(2)-EV_(2))+GU*(EV_(4)-2.0*EV_(3));
        GW = AW*(pbm.efpar(2)-EV_(2))-A+pbm.efpar(1)*(EV_(4)-2.0*EV_(3));
        GP = DD1*(pbm.efpar(2)-EV_(2))-2.0*GU;
        GPP = DD2*(pbm.efpar(2)-EV_(2));
        GWW = GPP-2.0*AW;
        GWP = GPP-DD1-2.0*pbm.efpar(1);
        varargout{1} = EV_(1)*F*G;
        if(nargout>1)
            g_(1,1) = F*G;
            g_(2,1) = EV_(1)*(FW*G+F*GW);
            g_(3,1) = EV_(1)*(FP*G+F*GP);
            g_(4,1) = EV_(1)*F*GU;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = FW*G+F*GW;
                H_(2,1) = H_(1,2);
                H_(1,3) = FP*G+F*GP;
                H_(3,1) = H_(1,3);
                H_(1,4) = F*GU;
                H_(4,1) = H_(1,4);
                H_(2,2) = EV_(1)*(FWW*G+2.0*FW*GW+F*GWW);
                H_(2,3) = EV_(1)*(FWP*G+FW*GP+FP*GW+F*GWP);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*(FW*GU+F*pbm.efpar(1));
                H_(4,2) = H_(2,4);
                H_(3,3) = EV_(1)*(FPP*G+2.0*FP*GP+F*GPP);
                H_(3,4) = EV_(1)*FP*GU;
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    case 'eEP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        D = 1.0+pbm.efpar(3)*(EV_(2)+EV_(3));
        DD = D*D;
        DD1 = 2.0*pbm.efpar(3)*D;
        DD2 = 2.0*pbm.efpar(3)*pbm.efpar(3);
        A = DD+pbm.efpar(4)+pbm.efpar(1)*EV_(2);
        AW = DD1+pbm.efpar(1);
        B = DD+pbm.efpar(4)+pbm.efpar(1)*EV_(3);
        BP = DD1+pbm.efpar(1);
        C = A*B-pbm.efpar(5)*EV_(2)*EV_(3);
        CW = AW*B+A*DD1-pbm.efpar(5)*EV_(3);
        CP = DD1*B+A*BP-pbm.efpar(5)*EV_(2);
        CWW = DD2*B+2.0*AW*DD1+A*DD2;
        CWP = DD2*B+AW*BP+DD1*DD1+A*DD2-pbm.efpar(5);
        CPP = DD2*B+2.0*DD1*BP+A*DD2;
        F = DD/C;
        H = DD1-F*CW;
        I = DD1-F*CP;
        FW = H/C;
        FP = I/C;
        HW = DD2-CW*FW-F*CWW;
        HP = DD2-CW*FW-F*CWP;
        IP = DD2-CP*FP-F*CPP;
        FWW = (HW-FW*CW)/C;
        FWP = (HP-FW*CP)/C;
        FPP = (IP-FP*CP)/C;
        G = B*(EV_(4)-2.0*EV_(3))+pbm.efpar(1)*EV_(3)*(pbm.efpar(2)-EV_(2));
        GW = DD1*(EV_(4)-2.0*EV_(3))-pbm.efpar(1)*EV_(3);
        GP = BP*(EV_(4)-2.0*EV_(3))-2.0*B+pbm.efpar(1)*(pbm.efpar(2)-EV_(2));
        GWW = DD2*(EV_(4)-2.0*EV_(3));
        GWP = GWW-2.0*DD1-pbm.efpar(1);
        GPP = GWW-4.0*BP;
        varargout{1} = EV_(1)*F*G;
        if(nargout>1)
            g_(1,1) = F*G;
            g_(2,1) = EV_(1)*(FW*G+F*GW);
            g_(3,1) = EV_(1)*(FP*G+F*GP);
            g_(4,1) = EV_(1)*F*B;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = FW*G+F*GW;
                H_(2,1) = H_(1,2);
                H_(1,3) = FP*G+F*GP;
                H_(3,1) = H_(1,3);
                H_(1,4) = F*B;
                H_(4,1) = H_(1,4);
                H_(2,2) = EV_(1)*(FWW*G+2.0*FW*GW+F*GWW);
                H_(2,3) = EV_(1)*(FWP*G+FW*GP+FP*GW+F*GWP);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*(FW*B+F*DD1);
                H_(4,2) = H_(2,4);
                H_(3,3) = EV_(1)*(FPP*G+2.0*FP*GP+F*GPP);
                H_(3,4) = EV_(1)*(FP*B+F*BP);
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [5,0];
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

