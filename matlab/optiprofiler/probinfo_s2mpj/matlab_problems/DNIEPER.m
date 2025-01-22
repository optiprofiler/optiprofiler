function varargout = DNIEPER(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DNIEPER
%    *********
% 
%    This problem models the planning of systematic use of water resources
%    in the basin of the river Dnieper.
% 
%    Source: p. 139sq in 
%    B.N. Pshenichnyj
%    "The Linearization Method for Constrained Optimization",
%    Springer Verlag, SCM Series 22, Heidelberg, 1994
% 
%    SIF input: Ph. Toint, December 1994.
% 
%    classification = 'C-CQOR2-MN-61-24'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DNIEPER';

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
        v_('C1') = 5.61;
        v_('C2') = 4.68;
        v_('C3') = 1.62;
        v_('C4') = 1.8;
        v_('C5') = 2.13;
        v_('C6') = 2.1;
        v_('C7') = 1.99;
        v_('C8') = 2.02;
        v_('C9') = 2.14;
        v_('C10') = 2.15;
        v_('C11') = 2.36;
        v_('C12') = 2.63;
        v_('C13') = -0.02;
        v_('C14') = -0.01;
        v_('C15') = -0.16;
        v_('C16') = -0.47;
        v_('C17') = -0.75;
        v_('C18') = -0.94;
        v_('C19') = -0.93;
        v_('C20') = -0.99;
        v_('C21') = -0.42;
        v_('C22') = -0.07;
        v_('C23') = 0.04;
        v_('C24') = -0.06;
        v_('1') = 1;
        v_('2') = 2;
        v_('4') = 4;
        v_('5') = 5;
        v_('8') = 8;
        v_('9') = 9;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_('16') = 16;
        v_('17') = 17;
        v_('20') = 20;
        v_('21') = 21;
        v_('24') = 24;
        v_('25') = 25;
        v_('36') = 36;
        v_('37') = 37;
        v_('48') = 48;
        v_('49') = 49;
        v_('52') = 52;
        v_('53') = 53;
        v_('56') = 56;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('56')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','X0F',ix_);
        pb.xnames{iv} = 'X0F';
        [iv,ix_] = s2mpjlib('ii','X24F',ix_);
        pb.xnames{iv} = 'X24F';
        [iv,ix_] = s2mpjlib('ii','X12F',ix_);
        pb.xnames{iv} = 'X12F';
        [iv,ix_] = s2mpjlib('ii','X36F',ix_);
        pb.xnames{iv} = 'X36F';
        [iv,ix_] = s2mpjlib('ii','AC',ix_);
        pb.xnames{iv} = 'AC';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('12')
            v_('I0') = 12+I;
            v_('I1') = 24+I;
            v_('I2') = 36+I;
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I1')))]);
            valA(end+1) = 19.95;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 0.07656;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I2')))]);
            valA(end+1) = -24.89;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -0.7135;
        end
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        pbm.gscale(ig,1) = -1.0;
        for I=v_('1'):v_('4')
            v_('I0') = 24+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
        end
        for I=v_('5'):v_('8')
            v_('I0') = 24+I;
            v_('I1') = 44+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I1')))]);
            valA(end+1) = -2.68;
        end
        for I=v_('9'):v_('12')
            v_('I0') = 24+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
        end
        for I=v_('13'):v_('16')
            v_('I0') = 12+I;
            v_('I1') = 24+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I1')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('AC');
            valA(end+1) = -1.0;
        end
        for I=v_('17'):v_('20')
            v_('I0') = 12+I;
            v_('I1') = 24+I;
            v_('I2') = 36+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I1')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I2')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('AC');
            valA(end+1) = -1.0;
        end
        for I=v_('21'):v_('24')
            v_('I0') = 12+I;
            v_('I1') = 24+I;
            [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CC',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I0')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I1')))]);
            valA(end+1) = -2.68;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('AC');
            valA(end+1) = -1.0;
        end
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
        pbm.gconst(ig_('OBJ')) = -112.464;
        for I=v_('1'):v_('24')
            v_('CST') = -1.0*v_(['C',int2str(I)]);
            pbm.gconst(ig_(['CC',int2str(I)])) = v_('CST');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('12')
            pb.xlower(ix_(['X',int2str(I)]),1) = 51.2;
            pb.xupper(ix_(['X',int2str(I)])) = 51.4;
        end
        for I=v_('13'):v_('24')
            pb.xlower(ix_(['X',int2str(I)]),1) = 15.0;
            pb.xupper(ix_(['X',int2str(I)])) = 16.1;
        end
        for I=v_('25'):v_('36')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.4;
            pb.xupper(ix_(['X',int2str(I)])) = 4.6;
        end
        for I=v_('37'):v_('48')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.5;
            pb.xupper(ix_(['X',int2str(I)])) = 4.8;
        end
        for I=v_('49'):v_('56')
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I)])) = 0.7;
        end
        pb.xlower(ix_('X0F'),1) = 50.82;
        pb.xupper(ix_('X0F'),1) = 50.82;
        pb.xlower(ix_('X24F'),1) = 2.0;
        pb.xupper(ix_('X24F'),1) = 2.0;
        pb.xlower(ix_('X12F'),1) = 15.5;
        pb.xupper(ix_('X12F'),1) = 15.5;
        pb.xlower(ix_('X36F'),1) = 2.3;
        pb.xupper(ix_('X36F'),1) = 2.3;
        pb.xlower(ix_('AC')) = -Inf;
        pb.xupper(ix_('AC'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('12')
            pb.x0(ix_(['X',int2str(I)]),1) = 51.35;
        end
        for I=v_('13'):v_('24')
            pb.x0(ix_(['X',int2str(I)]),1) = 15.5;
        end
        for I=v_('25'):v_('36')
            pb.x0(ix_(['X',int2str(I)]),1) = 2.5;
        end
        for I=v_('37'):v_('48')
            pb.x0(ix_(['X',int2str(I)]),1) = 2.6;
        end
        for I=v_('49'):v_('56')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.3;
        end
        pb.x0(ix_('X0F'),1) = 50.82;
        pb.x0(ix_('X24F'),1) = 2.0;
        pb.x0(ix_('X12F'),1) = 15.5;
        pb.x0(ix_('X36F'),1) = 2.3;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eWJ',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eWK',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('12')
            v_('I1') = 12+I;
            v_('I2') = 36+I;
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'ACSQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQ';
        ielftype(ie) = iet_('eSQ');
        vname = 'AC';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('1'):v_('12')
            v_('I0') = 24+I;
            ename = ['W1',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWJ';
            ielftype(ie) = iet_('eWJ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['W2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eWJ';
        ielftype(ie) = iet_('eWJ');
        ename = ['W2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'X0F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'X24F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('2'):v_('12')
            v_('I0') = 23+I;
            v_('I1') = -1+I;
            ename = ['W2',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWJ';
            ielftype(ie) = iet_('eWJ');
            vname = ['X',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('13'):v_('24')
            v_('I1') = 24+I;
            ename = ['W1',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWK';
            ielftype(ie) = iet_('eWK');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = ['W2',int2str(round(v_('13')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eWK';
        ielftype(ie) = iet_('eWK');
        ename = ['W2',int2str(round(v_('13')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'X12F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['W2',int2str(round(v_('13')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'X36F';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('14'):v_('24')
            v_('I0') = 23+I;
            v_('I1') = -1+I;
            ename = ['W2',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eWK';
            ielftype(ie) = iet_('eWK');
            vname = ['X',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I0')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('12')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 2.155;
        end
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('ACSQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2000.0;
        for I=v_('1'):v_('24')
            ig = ig_(['CC',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['W1',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['W2',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLUTION             1.87439D+04
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MN-61-24';
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

    case 'eWJ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A1 = 34.547;
        A2 = -0.55878;
        A3 = 8.05339;
        A4 = -0.02252;
        A5 = -0.29316;
        A6 = -0.013521;
        A7 = 0.00042;
        A8 = 0.00267;
        A9 = 0.000281;
        A10 = 0.0000032;
        varargout{1} = A1+A2*EV_(1)+A3*EV_(2)+A4*EV_(1)^2+A5*EV_(1)*EV_(2)+...
             A6*EV_(2)^2+A7*EV_(1)^3+A8*EV_(1)^2*EV_(2)+A9*EV_(1)*EV_(2)^2+A10*EV_(2)^3;
        if(nargout>1)
            g_(1,1) = A2+2.0*A4*EV_(1)+A5*EV_(2)+3.0*A7*EV_(1)^2+2.0*A8*EV_(1)*...
                 EV_(2)+A9*EV_(2)^2;
            g_(2,1) = A3+A5*EV_(1)+2.0*A6*EV_(2)+A8*EV_(1)^2+2.0*A9*EV_(1)*EV_(2)+...
                 3.0*A10*EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*A4+6.0*A7*EV_(1)+2.0*A8*EV_(2);
                H_(1,2) = A5+2.0*A8*EV_(1)+2.0*A9*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*A6+2.0*A9*EV_(1)+6.0*A10*EV_(2);
                varargout{3} = H_;
            end
        end

    case 'eWK'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A1 = 20.923;
        A2 = -4.22088;
        A3 = 1.42061;
        A4 = -0.41040;
        A5 = -0.15082;
        A7 = -0.00826;
        A8 = 0.00404;
        A9 = 0.000168;
        A10 = -0.000038;
        varargout{1} = A1+A2*EV_(1)+A3*EV_(2)+A4*EV_(1)^2+A5*EV_(1)*EV_(2)+...
             A7*EV_(1)^3+A8*EV_(1)^2*EV_(2)+A9*EV_(1)*EV_(2)^2+A10*EV_(2)^3;
        if(nargout>1)
            g_(1,1) = A2+2.0*A4*EV_(1)+A5*EV_(2)+3.0*A7*EV_(1)^2+2.0*A8*EV_(1)*...
                 EV_(2)+A9*EV_(2)^2;
            g_(2,1) = A3+A5*EV_(1)+A8*EV_(1)^2+2.0*A9*EV_(1)*EV_(2)+3.0*A10*EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*A4+6.0*A7*EV_(1)+2.0*A8*EV_(2);
                H_(1,2) = A5+2.0*A8*EV_(1)+2.0*A9*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*A9*EV_(1)+6.0*A10*EV_(2);
                varargout{3} = H_;
            end
        end

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
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
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

