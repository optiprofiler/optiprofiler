function varargout = RES(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : RES
%    *********
% 
%    Dassault France ressort (spring) problem
% 
%    SIF input:  A. R. Conn, June 1993.
% 
%    classification = 'C-CNLR2-MN-20-14'
% 
% 
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'RES';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','L0',ix_);
        pb.xnames{iv} = 'L0';
        [iv,ix_] = s2mpjlib('ii','N',ix_);
        pb.xnames{iv} = 'N';
        [iv,ix_] = s2mpjlib('ii','F',ix_);
        pb.xnames{iv} = 'F';
        [iv,ix_] = s2mpjlib('ii','K',ix_);
        pb.xnames{iv} = 'K';
        [iv,ix_] = s2mpjlib('ii','LB',ix_);
        pb.xnames{iv} = 'LB';
        [iv,ix_] = s2mpjlib('ii','L',ix_);
        pb.xnames{iv} = 'L';
        [iv,ix_] = s2mpjlib('ii','DE',ix_);
        pb.xnames{iv} = 'DE';
        [iv,ix_] = s2mpjlib('ii','DI',ix_);
        pb.xnames{iv} = 'DI';
        [iv,ix_] = s2mpjlib('ii','TO',ix_);
        pb.xnames{iv} = 'TO';
        [iv,ix_] = s2mpjlib('ii','TOB',ix_);
        pb.xnames{iv} = 'TOB';
        [iv,ix_] = s2mpjlib('ii','NU',ix_);
        pb.xnames{iv} = 'NU';
        [iv,ix_] = s2mpjlib('ii','D',ix_);
        pb.xnames{iv} = 'D';
        [iv,ix_] = s2mpjlib('ii','P',ix_);
        pb.xnames{iv} = 'P';
        [iv,ix_] = s2mpjlib('ii','E',ix_);
        pb.xnames{iv} = 'E';
        [iv,ix_] = s2mpjlib('ii','P0',ix_);
        pb.xnames{iv} = 'P0';
        [iv,ix_] = s2mpjlib('ii','G',ix_);
        pb.xnames{iv} = 'G';
        [iv,ix_] = s2mpjlib('ii','DM',ix_);
        pb.xnames{iv} = 'DM';
        [iv,ix_] = s2mpjlib('ii','FR',ix_);
        pb.xnames{iv} = 'FR';
        [iv,ix_] = s2mpjlib('ii','TOLIM',ix_);
        pb.xnames{iv} = 'TOLIM';
        [iv,ix_] = s2mpjlib('ii','TOBLIM',ix_);
        pb.xnames{iv} = 'TOBLIM';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','E1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('K');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DE');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('D');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DM');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DI');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('D');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('DM');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('D');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('P');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('E');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('NU');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('N');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('D');
        valA(end+1) = 1.5;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('L0');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('L');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LB');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('FR');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E9';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('LB');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E10';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('L');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('L0');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TO');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'E12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TOB');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','E13',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'E13';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TO');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TOLIM');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','E14',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'E14';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TOB');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('TOBLIM');
        valA(end+1) = -1.0;
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
        pbm.gconst(ig_('E6')) = -2.0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('L0')) = 100.0;
        pb.xupper(ix_('N')) = 100.0;
        pb.xupper(ix_('F')) = 30.0;
        pb.xupper(ix_('K')) = 100.0;
        pb.xupper(ix_('LB')) = 50.0;
        pb.xupper(ix_('L')) = 50.0;
        pb.xupper(ix_('DE')) = 30.0;
        pb.xupper(ix_('DI')) = 30.0;
        pb.xupper(ix_('TO')) = 800.0;
        pb.xupper(ix_('TOB')) = 800.0;
        pb.xupper(ix_('NU')) = 50.0;
        pb.xlower(ix_('NU'),1) = 0.5;
        pb.xupper(ix_('D')) = 10.0;
        pb.xlower(ix_('D'),1) = 0.1;
        pb.xupper(ix_('P')) = 20.0;
        pb.xupper(ix_('E')) = 10.0;
        pb.xupper(ix_('P0')) = 1000.0;
        pb.xlower(ix_('P0'),1) = 1.0;
        pb.xupper(ix_('G')) = 80000.0;
        pb.xlower(ix_('G'),1) = 40000.0;
        pb.xupper(ix_('DM')) = 30.0;
        pb.xlower(ix_('DM'),1) = 0.1;
        pb.xupper(ix_('FR')) = 50.0;
        pb.xupper(ix_('TOLIM')) = 1000.0;
        pb.xlower(ix_('TOLIM'),1) = 100.0;
        pb.xupper(ix_('TOBLIM')) = 1000.0;
        pb.xlower(ix_('TOBLIM'),1) = 100.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('L0'),1) = 1.5000e-01;
        pb.x0(ix_('N'),1) = 2.4079e+01;
        pb.x0(ix_('F'),1) = 9.2459e-15;
        pb.x0(ix_('K'),1) = 0.0000;
        pb.x0(ix_('LB'),1) = 0.0000;
        pb.x0(ix_('L'),1) = 1.5000e-01;
        pb.x0(ix_('DE'),1) = 6.8120;
        pb.x0(ix_('DI'),1) = 6.6120;
        pb.x0(ix_('TO'),1) = 0.0000;
        pb.x0(ix_('TOB'),1) = 0.0000;
        pb.x0(ix_('NU'),1) = 2.2079e+01;
        pb.x0(ix_('D'),1) = 1.0000e-01;
        pb.x0(ix_('P'),1) = 6.5268e-01;
        pb.x0(ix_('E'),1) = 5.5268e-01;
        pb.x0(ix_('P0'),1) = 6.5887e+02;
        pb.x0(ix_('G'),1) = 6.5887e+04;
        pb.x0(ix_('DM'),1) = 6.7120;
        pb.x0(ix_('FR'),1) = 1.5000e-01;
        pb.x0(ix_('TOLIM'),1) = 1.0000e+02;
        pb.x0(ix_('TOBLIM'),1) = 1.0000e+02;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'en311d14',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        elftv{it}{3} = 'X';
        elftv{it}{4} = 'Y';
        elftv{it}{5} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en14d31',iet_);
        elftv{it}{1} = 'W';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        elftv{it}{4} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en11d3',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en111d2',iet_);
        elftv{it}{1} = 'W';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        elftv{it}{4} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'EL1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en311d14';
        ielftype(ie) = iet_('en311d14');
        vname = 'DM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NU';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'P0';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'G';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en14d31';
        ielftype(ie) = iet_('en14d31');
        vname = 'G';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'DM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NU';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'NU';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'P';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en11d3';
        ielftype(ie) = iet_('en11d3');
        vname = 'P0';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'DM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EL5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en111d2';
        ielftype(ie) = iet_('en111d2');
        vname = 'G';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'D';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'E';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'DM';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CNLR2-MN-20-14';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 3.1415926535;
        varargout{1} = pbm;

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

    case 'en311d14'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V3WX = EV_(1)^3*EV_(2)*EV_(3);
        YZ4 = EV_(4)*EV_(5)^4;
        varargout{1} = V3WX/YZ4;
        if(nargout>1)
            g_(1,1) = (3.0*EV_(1)^2*EV_(2)*EV_(3))/YZ4;
            g_(2,1) = (EV_(1)^3*EV_(3))/YZ4;
            g_(3,1) = (EV_(1)^3*EV_(2))/YZ4;
            g_(4,1) = -V3WX/(EV_(4)*YZ4);
            g_(5,1) = -(4.0*V3WX)/(EV_(5)*YZ4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = (6.0*EV_(1)*EV_(2)*EV_(3))/YZ4;
                H_(1,2) = (3.0*EV_(1)^2*EV_(3))/YZ4;
                H_(2,1) = H_(1,2);
                H_(1,3) = (3.0*EV_(1)^2*EV_(2))/YZ4;
                H_(3,1) = H_(1,3);
                H_(1,4) = -(3.0*EV_(1)^2*EV_(2))/(YZ4*EV_(4));
                H_(4,1) = H_(1,4);
                H_(1,5) = -(12.0*EV_(1)^2*EV_(2)*EV_(3))/(YZ4*EV_(5));
                H_(5,1) = H_(1,5);
                H_(2,3) = EV_(1)^3/YZ4;
                H_(3,2) = H_(2,3);
                H_(2,4) = -(EV_(1)^3*EV_(3))/(YZ4*EV_(4));
                H_(4,2) = H_(2,4);
                H_(2,5) = -(4.0*EV_(1)^3*EV_(3))/(YZ4*EV_(5));
                H_(5,2) = H_(2,5);
                H_(3,4) = -(EV_(1)^3*EV_(2))/(YZ4*EV_(4));
                H_(4,3) = H_(3,4);
                H_(3,5) = -(4.0*EV_(1)^3*EV_(2))/(YZ4*EV_(5));
                H_(5,3) = H_(3,5);
                H_(4,4) = -(2.0*V3WX)/(EV_(4)^2*YZ4);
                H_(4,5) = (4.0*V3WX)/(EV_(4)*EV_(5)*YZ4);
                H_(5,4) = H_(4,5);
                H_(5,5) = (20.0*V3WX)/(EV_(5)^2*YZ4);
                varargout{3} = H_;
            end
        end

    case 'en14d31'

        EV_  = varargin{1};
        iel_ = varargin{2};
        WX4 = EV_(1)*EV_(2)^4;
        Y3Z = EV_(3)^3*EV_(4);
        varargout{1} = WX4/Y3Z;
        if(nargout>1)
            g_(1,1) = EV_(2)^4/Y3Z;
            g_(2,1) = (4.0*EV_(1)*EV_(2)^3)/Y3Z;
            g_(3,1) = -(3.0*WX4)/(EV_(3)*Y3Z);
            g_(4,1) = -WX4/(EV_(4)*Y3Z);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = (4.0*EV_(2)^3)/Y3Z;
                H_(2,1) = H_(1,2);
                H_(1,3) = -(3.0*EV_(2)^4)/(EV_(3)*Y3Z);
                H_(3,1) = H_(1,3);
                H_(1,4) = -EV_(2)^4/(EV_(4)*Y3Z);
                H_(4,1) = H_(1,4);
                H_(2,2) = (12.0*EV_(1)*EV_(2)^2)/Y3Z;
                H_(2,3) = -(12.0*EV_(1)*EV_(2)^3)/(EV_(3)*Y3Z);
                H_(3,2) = H_(2,3);
                H_(2,4) = -(4.0*EV_(1)*EV_(2)^3)/(EV_(4)*Y3Z);
                H_(4,2) = H_(2,4);
                H_(3,3) = (12.0*WX4)/(EV_(3)^2*Y3Z);
                H_(3,4) = (3.0*WX4)/(EV_(4)*EV_(3)*Y3Z);
                H_(4,3) = H_(3,4);
                H_(4,4) = (2.0*WX4)/(EV_(4)^2*Y3Z);
                varargout{3} = H_;
            end
        end

    case 'en11d3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)*EV_(2))/(pbm.efpar(1)*EV_(3)^3);
        if(nargout>1)
            g_(1,1) = EV_(2)/(pbm.efpar(1)*EV_(3)^3);
            g_(2,1) = EV_(1)/(pbm.efpar(1)*EV_(3)^3);
            g_(3,1) = -(3.0*EV_(1)*EV_(2))/(pbm.efpar(1)*EV_(3)^4);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = 1.0/(pbm.efpar(1)*EV_(3)^3);
                H_(2,1) = H_(1,2);
                H_(1,3) = -(3.0*EV_(2))/(pbm.efpar(1)*EV_(3)^4);
                H_(3,1) = H_(1,3);
                H_(2,3) = -(3.0*EV_(1))/(pbm.efpar(1)*EV_(3)^4);
                H_(3,2) = H_(2,3);
                H_(3,3) = (12.0*EV_(1)*EV_(2))/(pbm.efpar(1)*EV_(3)^5);
                varargout{3} = H_;
            end
        end

    case 'en111d2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)*EV_(2)*EV_(3))/(pbm.efpar(1)*EV_(4)^2);
        if(nargout>1)
            g_(1,1) = (EV_(2)*EV_(3))/(pbm.efpar(1)*EV_(4)^2);
            g_(2,1) = (EV_(1)*EV_(3))/(pbm.efpar(1)*EV_(4)^2);
            g_(3,1) = (EV_(1)*EV_(2))/(pbm.efpar(1)*EV_(4)^2);
            g_(4,1) = -(2.0*EV_(1)*EV_(2)*EV_(3))/(pbm.efpar(1)*EV_(4)^3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = EV_(3)/(pbm.efpar(1)*EV_(4)^2);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)/(pbm.efpar(1)*EV_(4)^2);
                H_(3,1) = H_(1,3);
                H_(1,4) = -(2.0*EV_(2)*EV_(3))/(pbm.efpar(1)*EV_(4)^3);
                H_(4,1) = H_(1,4);
                H_(2,3) = EV_(1)/(pbm.efpar(1)*EV_(4)^2);
                H_(3,2) = H_(2,3);
                H_(2,4) = -(2.0*EV_(1)*EV_(3))/(pbm.efpar(1)*EV_(4)^3);
                H_(4,2) = H_(2,4);
                H_(3,4) = -(2.0*EV_(1)*EV_(2))/(pbm.efpar(1)*EV_(4)^3);
                H_(4,3) = H_(3,4);
                H_(4,4) = (6.0*EV_(1)*EV_(2)*EV_(3))/(pbm.efpar(1)*EV_(4)^4);
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

