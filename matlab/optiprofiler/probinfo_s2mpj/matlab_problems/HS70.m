function varargout = HS70(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS70
%    *********
% 
%    This problem arises in water flow routing.
% 
%    Source: problem 70 incorrectly stated in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991, modified May 2024
% 
%    classification = 'C-CSQR2-MN-4-1'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS70';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 4;
        v_('1') = 1;
        v_('2') = 2;
        v_('19') = 19;
        v_('C1') = 0.1;
        v_('C2') = 1.0;
        v_('C3') = 2.0;
        v_('C4') = 3.0;
        v_('C5') = 4.0;
        v_('C6') = 5.0;
        v_('C7') = 6.0;
        v_('C8') = 7.0;
        v_('C9') = 8.0;
        v_('C10') = 9.0;
        v_('C11') = 10.0;
        v_('C12') = 11.0;
        v_('C13') = 12.0;
        v_('C14') = 13.0;
        v_('C15') = 14.0;
        v_('C16') = 15.0;
        v_('C17') = 16.0;
        v_('C18') = 17.0;
        v_('C19') = 18.0;
        v_('Y1') = 0.00189;
        v_('Y2') = 0.1038;
        v_('Y3') = 0.268;
        v_('Y4') = 0.506;
        v_('Y5') = 0.577;
        v_('Y6') = 0.604;
        v_('Y7') = 0.725;
        v_('Y8') = 0.898;
        v_('Y9') = 0.947;
        v_('Y10') = 0.845;
        v_('Y11') = 0.702;
        v_('Y12') = 0.528;
        v_('Y13') = 0.385;
        v_('Y14') = 0.257;
        v_('Y15') = 0.159;
        v_('Y16') = 0.0869;
        v_('Y17') = 0.0453;
        v_('Y18') = 0.01509;
        v_('Y19') = 0.00189;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('19')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 1.0e+0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 1.0e+0;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        for I=v_('1'):v_('19')
            pbm.gconst(ig_(['OBJ',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.00001*ones(pb.n,1);
        pb.xupper = 100.0*ones(pb.n,1);
        pb.xupper(ix_('X3')) = 1.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 2.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 4.0;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 4.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.04;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.04;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eY1',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eY2',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X3';
        elftv{it}{2} = 'X4';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('19')
            ename = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eY1';
            ielftype(ie) = iet_('eY1');
            ename = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('C',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(I)]);
            ename = ['Y',int2str(I),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eY2';
            ielftype(ie) = iet_('eY2');
            ename = ['Y',int2str(I),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['Y',int2str(I),',',int2str(round(v_('2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            [~,posep] = ismember('C',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(I)]);
        end
        ename = 'C1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.00001,100.0,[]);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gSQR',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('19')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gSQR';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Y',int2str(I),',',int2str(round(v_('1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Y',int2str(I),',',int2str(round(v_('2')))]);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.007498464
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CSQR2-MN-4-1';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
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
        pbm.efpar(1) = sqrt(1.0e+0/6.2832e+0);
        varargout{1} = pbm;

    case 'eY1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(2)+EV_(3)*(1.0e+0-EV_(2));
        CI = pbm.elpar{iel_}(1)/7.658e+0;
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_(1));
        P0V1 = -1.0e+0/(1.2e+1*EV_(1)^2);
        P0V1V1 = 2.0e+0/(1.2e+1*EV_(1)^3);
        P1 = 1.0e+0/P0;
        P1V1 = -P0V1/(P0^2);
        P1V1V1 = (2.0e+0*P0V1^2/P0-P0V1V1)/(P0^2);
        P2 = EV_(2);
        P3 = B^EV_(1);
        LOGB = log(B);
        P3V1 = P3*LOGB;
        P3V2 = EV_(1)*(1.0e+0-EV_(3))*B^(EV_(1)-1.0e+0);
        P3V3 = EV_(1)*(1.0e+0-EV_(2))*B^(EV_(1)-1.0e+0);
        P3V1V1 = P3V1*LOGB;
        P3V1V2 = P3V2*LOGB+P3*(1.0e+0-EV_(3))/B;
        P3V1V3 = P3V3*LOGB+P3*(1.0e+0-EV_(2))/B;
        P3V2V2 = EV_(1)*(EV_(1)-1.0e+0)*(1.0e+0-EV_(3))^2*B^(EV_(1)-1.0e+0);
        P3V2V3 = -EV_(1)*B^(EV_(1)-1.0e+0)+EV_(1)*(EV_(1)-1.0e+0)*(1.0e+0-EV_(2))*...
             (1.0e+0-EV_(3))*B^(EV_(1)-2.0e+0);
        P3V3V3 = EV_(1)*(EV_(1)-1.0e+0)*(1.0e+0-EV_(2))^2*B^(EV_(1)-2.0e+0);
        P4 = pbm.efpar(1)*sqrt(EV_(1));
        P4V1 = 5.0e-1*pbm.efpar(1)*sqrt(1.0e+0/EV_(1));
        P4V1V1 = -2.5e-1*pbm.efpar(1)*sqrt(1.0e+0/EV_(1)^3);
        C5 = CI^(-1.0e+0);
        P5 = C5*CI^EV_(1);
        P5V1 = P5*log(CI);
        P5V1V1 = P5V1*log(CI);
        P6 = exp(EV_(1)*(1.0e+0-CI*B));
        P6V1 = P6*(1.0e+0-CI*B);
        P6V2 = -P6*EV_(1)*CI*(1.0e+0-EV_(3));
        P6V3 = -P6*EV_(1)*CI*(1.0e+0-EV_(2));
        P6V1V1 = P6*(1.0e+0-CI*B)^2;
        P6V1V2 = P6V2*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_(3));
        P6V1V3 = P6V3*(1.0e+0-CI*B)-P6*CI*(1.0e+0-EV_(2));
        P6V2V2 = -P6V2*EV_(1)*CI*(1.0e+0-EV_(3));
        P6V2V3 = -P6V3*EV_(1)*CI*(1.0e+0-EV_(3))+P6*EV_(1)*CI;
        P6V3V3 = -P6V3*EV_(1)*CI*(1.0e+0-EV_(2));
        varargout{1} = P1*P2*P3*P4*P5*P6;
        if(nargout>1)
            g_(1,1) = P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*...
                 P6+P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1;
            g_(2,1) = P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2;
            g_(3,1) = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*...
                     P5*P6+P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1);
                H_(1,2) = P1V1*(P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*...
                     (P3V1*P4*P5*P6+(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2)));
                H_(2,1) = H_(1,2);
                H_(1,3) =...
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3)));
                H_(3,1) = H_(1,3);
                H_(2,2) =...
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2+P3*P6V2+P3V2*P6));
                H_(2,3) =...
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3+P3V3*P6+P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2);
                H_(3,2) = H_(2,3);
                H_(3,3) = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3);
                varargout{3} = H_;
            end
        end

    case 'eY2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        B = EV_(2)+EV_(3)*(1.0e+0-EV_(2));
        CI = pbm.elpar{iel_}(1)/7.658e+0;
        P0 = 1.0e+0+1.0e+0/(1.2e+1*EV_(1));
        P0V1 = -1.0e+0/(1.2e+1*EV_(1)^2);
        P0V1V1 = 2.0e+0/(1.2e+1*EV_(1)^3);
        P1 = 1.0e+0/P0;
        P1V1 = -P0V1/(P0^2);
        P1V1V1 = (2.0e+0*P0V1^2/P0-P0V1V1)/(P0^2);
        P2 = 1.0e+0-EV_(2);
        P3 = (B/EV_(3))^EV_(1);
        LOGB = log(B/EV_(3));
        P3V1 = P3*LOGB;
        P3V2 = EV_(1)*(-1.0e+0+1.0e+0/EV_(3))*(B/EV_(3))^(EV_(1)-1.0e+0);
        P3V3 = -EV_(1)*(EV_(2)/EV_(3)^2)*(B/EV_(3))^(EV_(1)-1.0e+0);
        P3V1V1 = P3V1*LOGB;
        P3V1V2 = P3V2*LOGB+P3*EV_(3)*(-1.0e+0+1.0e+0/EV_(3))/B;
        P3V1V3 = P3V3*LOGB-P3*EV_(2)/(B*EV_(3));
        P3V2V2 =...
              EV_(1)*(EV_(1)-1.0e+0)*(-1.0e+0+1.0e+0/EV_(3))^2*(B/EV_(3))^(EV_(1)-2.0e+0);
        P3V2V3 = EV_(1)*(-1.0e+0/EV_(3)^2)*(B/EV_(3))^(EV_(1)-1.0e+0)+EV_(1)*...
             (EV_(1)-1.0e+0)*(-1.0e+0+1.0e+0/EV_(3))*(-EV_(2)/EV_(3)^2)*(B/EV_(3))^(EV_(1)-2.0e+0);
        P3V3V3 = 2.0e+0*EV_(1)*(EV_(2)/EV_(3)^3)*(B/EV_(3))^(EV_(1)-1.0e+0)+...
             EV_(1)*(EV_(1)-1.0e+0)*(EV_(2)/EV_(3)^2)^2*(B/EV_(3))^(EV_(1)-2.0e+0);
        P4 = pbm.efpar(1)*sqrt(EV_(1));
        P4V1 = 5.0e-1*pbm.efpar(1)*sqrt(1.0e+0/EV_(1));
        P4V1V1 = -2.5e-1*pbm.efpar(1)*sqrt(1.0e+0/EV_(1)^3);
        C5 = CI^(-1.0e+0);
        P5 = C5*CI^EV_(1);
        P5V1 = P5*log(CI);
        P5V1V1 = P5V1*log(CI);
        P6 = exp(EV_(1)*(1.0e+0-CI*B/EV_(3)));
        P6V1 = P6*(1.0e+0-CI*B/EV_(3));
        P6V2 = -P6*EV_(1)*CI*(1.0e+0-EV_(3))/EV_(3);
        P6V3 = P6*EV_(1)*CI*EV_(2)/EV_(3)^2;
        P6V1V1 = P6*(1.0e+0-CI*B/EV_(3))^2;
        P6V1V2 = P6V2*(1.0e+0-CI*B/EV_(3))-P6*CI*(-1.0e+0+1.0e+0/EV_(3));
        P6V1V3 = P6V3*(1.0e+0-CI*B/EV_(3))+P6*CI*EV_(2)/EV_(3)^2;
        P6V2V2 = -P6V2*EV_(1)*CI*(1.0e+0-EV_(3))/EV_(3);
        P6V2V3 = -P6V3*EV_(1)*CI*(1.0e+0-EV_(3))/EV_(3)+P6*EV_(1)*CI/EV_(3)^2;
        P6V3V3 = P6V3*EV_(1)*CI*EV_(2)/EV_(3)^2-2.0e+0*P6*EV_(1)*CI*EV_(2)/EV_(3)^3;
        varargout{1} = P1*P2*P3*P4*P5*P6;
        if(nargout>1)
            g_(1,1) = P1V1*P2*P3*P4*P5*P6+P1*P2*P3V1*P4*P5*P6+P1*P2*P3*P4V1*P5*...
                 P6+P1*P2*P3*P4*P5V1*P6+P1*P2*P3*P4*P5*P6V1;
            g_(2,1) = -P1*P3*P4*P5*P6+P1*P2*P3V2*P4*P5*P6+P1*P2*P3*P4*P5*P6V2;
            g_(3,1) = P1*P2*P3V3*P4*P5*P6+P1*P2*P3*P4*P5*P6V3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = P1V1V1*P2*P3*P4*P5*P6+P1*P2*P3V1V1*P4*P5*P6+P1*P2*P3*P4V1V1*...
                     P5*P6+P1*P2*P3*P4*P5V1V1*P6+P1*P2*P3*P4*P5*P6V1V1+2.0e+0*(P1V1*P2*P3V1*P4*P5*P6+P1V1*P2*P3*P4V1*P5*P6+P1V1*P2*P3*P4*P5V1*P6+P1V1*P2*P3*P4*P5*P6V1+P1*P2*P3V1*P4V1*P5*P6+P1*P2*P3V1*P4*P5V1*P6+P1*P2*P3V1*P4*P5*P6V1+P1*P2*P3*P4V1*P5V1*P6+P1*P2*P3*P4V1*P5*P6V1+P1*P2*P3*P4*P5V1*P6V1);
                H_(1,2) = P1V1*(-P3*P4*P5*P6+P2*P3V2*P4*P5*P6+P2*P3*P4*P5*P6V2)+P1*...
                     (-P3V1*P4*P5*P6-(P3*P4V1*P5*P6+P3*P4*P5V1*P6+P3*P4*P5*P6V1)+P2*(P3V1V2*P4*P5*P6+P3V1*P4*P5*P6V2+P3V2*P4V1*P5*P6+P3V2*P4*P5V1*P6+P3V2*P4*P5*P6V1+P3*(P4V1*P5*P6V2+P4*P5V1*P6V2+P4*P5*P6V1V2)));
                H_(2,1) = H_(1,2);
                H_(1,3) =...
                      P2*(P1V1*P3V3*P4*P5*P6+P1V1*P3*P4*P5*P6V3+P1*(P3V1V3*P4*P5*P6+P3V1*P4*P5*P6V3+P3V3*P4V1*P5*P6+P3V3*P4*P5V1*P6+P3V3*P4*P5*P6V1+P3*(P4*P5V1*P6V3+P4V1*P5*P6V3+P4*P5*P6V1V3)));
                H_(3,1) = H_(1,3);
                H_(2,2) =...
                      P1*P4*P5*(P2*P3*P6V2V2+P2*P3V2V2*P6+2.0e+0*(P2*P3V2*P6V2-P3*P6V2-P3V2*P6));
                H_(2,3) =...
                      P1*P4*P5*(P2*P3V2V3*P6+P2*P3*P6V2V3-P3V3*P6-P3*P6V3+P2*P3V2*P6V3+P2*P3V3*P6V2);
                H_(3,2) = H_(2,3);
                H_(3,3) = P1*P2*P4*P5*(P3V3V3*P6+P3*P6V3V3+2.0e+0*P3V3*P6V3);
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
                H_(1,2) = 1.0e+0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQR'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0e+0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

