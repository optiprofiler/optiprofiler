function varargout = DEMBO7(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DEMBO7
%    *******
% 
%    A 7 stage membrane separation model
% 
%    Source: problem 7 in
%    R.S. Dembo,
%    "A set of geometric programming test problems and their solutions",
%    Mathematical Programming, 17, 192-213, 1976.
% 
%    SIF input: A. R. Conn, June 1993.
% 
%    classification = 'C-CQOR2-MN-16-20'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DEMBO7';

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
        v_('N') = 16;
        v_('1') = 1;
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
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 1.262626;
        [ig,ig_] = s2mpjlib('ii','C0',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'C0';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 1.262626;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 1.262626;
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 0.975;
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.975;
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.975;
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.975;
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 0.975;
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C6';
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = -0.002;
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C8';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','C9',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C9';
        [ig,ig_] = s2mpjlib('ii','C10',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C10';
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C11';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 0.002;
        [ig,ig_] = s2mpjlib('ii','C12',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C12';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 0.002;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = -0.002;
        [ig,ig_] = s2mpjlib('ii','C13',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C13';
        [ig,ig_] = s2mpjlib('ii','C14',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C14';
        [ig,ig_] = s2mpjlib('ii','C15',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C15';
        [ig,ig_] = s2mpjlib('ii','C16',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C16';
        [ig,ig_] = s2mpjlib('ii','C17',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C17';
        [ig,ig_] = s2mpjlib('ii','C18',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C18';
        [ig,ig_] = s2mpjlib('ii','C19',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'C19';
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
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0*ones(ngrp,1);
        pbm.gconst(ig_('OBJ')) = 0.0;
        pbm.gconst(ig_('C0')) = 50.0;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        grange(ig_('C0')) = 200.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.1*ones(pb.n,1);
        pb.xupper = 0.9*ones(pb.n,1);
        pb.xupper(ix_('X5')) = 1.0;
        pb.xlower(ix_('X5'),1) = 0.9;
        pb.xupper(ix_('X6')) = 0.1;
        pb.xlower(ix_('X6'),1) = 0.0001;
        pb.xupper(ix_('X11')) = 1000.0;
        pb.xlower(ix_('X11'),1) = 1.0;
        pb.xupper(ix_('X12')) = 500.0;
        pb.xlower(ix_('X12'),1) = 0.000001;
        pb.xupper(ix_('X13')) = 500.0;
        pb.xlower(ix_('X13'),1) = 1.0;
        pb.xupper(ix_('X14')) = 1000.0;
        pb.xlower(ix_('X14'),1) = 500.0;
        pb.xupper(ix_('X15')) = 1000.0;
        pb.xlower(ix_('X15'),1) = 500.0;
        pb.xupper(ix_('X16')) = 500.0;
        pb.xlower(ix_('X16'),1) = 0.000001;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 0.8;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 0.8;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 0.83;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 0.83;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 0.85;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 0.85;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 0.87;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 0.87;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 0.90;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 0.90;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 0.10;
        else
            pb.y0(find(pbm.congrps==ig('X6')),1) = 0.10;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 0.12;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 0.12;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 0.19;
        else
            pb.y0(find(pbm.congrps==ig('X8')),1) = 0.19;
        end
        if(isKey(ix_,'X9'))
            pb.x0(ix_('X9'),1) = 0.25;
        else
            pb.y0(find(pbm.congrps==ig_('X9')),1) = 0.25;
        end
        if(isKey(ix_,'X10'))
            pb.x0(ix_('X10'),1) = 0.29;
        else
            pb.y0(find(pbm.congrps==ig('X10')),1) = 0.29;
        end
        if(isKey(ix_,'X11'))
            pb.x0(ix_('X11'),1) = 512.0;
        else
            pb.y0(find(pbm.congrps==ig_('X11')),1) = 512.0;
        end
        if(isKey(ix_,'X12'))
            pb.x0(ix_('X12'),1) = 13.1;
        else
            pb.y0(find(pbm.congrps==ig('X12')),1) = 13.1;
        end
        if(isKey(ix_,'X13'))
            pb.x0(ix_('X13'),1) = 71.8;
        else
            pb.y0(find(pbm.congrps==ig_('X13')),1) = 71.8;
        end
        if(isKey(ix_,'X14'))
            pb.x0(ix_('X14'),1) = 640.0;
        else
            pb.y0(find(pbm.congrps==ig('X14')),1) = 640.0;
        end
        if(isKey(ix_,'X15'))
            pb.x0(ix_('X15'),1) = 650.0;
        else
            pb.y0(find(pbm.congrps==ig_('X15')),1) = 650.0;
        end
        if(isKey(ix_,'X16'))
            pb.x0(ix_('X16'),1) = 5.7;
        else
            pb.y0(find(pbm.congrps==ig('X16')),1) = 5.7;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eINV',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eQT',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eQTQT',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'W';
        elftv{it}{4} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'en2PRRC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eQTRC',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftv{it}{3} = 'Z';
        [it,iet_] = s2mpjlib( 'ii', 'eSQQT',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQT';
        ielftype(ie) = iet_('eSQQT');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQT';
        ielftype(ie) = iet_('eSQQT');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQT';
        ielftype(ie) = iet_('eSQQT');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQT';
        ielftype(ie) = iet_('eSQQT');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E14';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E15';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQQT';
        ielftype(ie) = iet_('eSQQT');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E16';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTQT';
        ielftype(ie) = iet_('eQTQT');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTQT';
        ielftype(ie) = iet_('eQTQT');
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PRRC';
        ielftype(ie) = iet_('en2PRRC');
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E21';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PRRC';
        ielftype(ie) = iet_('en2PRRC');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E22';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PRRC';
        ielftype(ie) = iet_('en2PRRC');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E23';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E24';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E25';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E26';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PR';
        ielftype(ie) = iet_('en2PR');
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E27';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E28';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTQT';
        ielftype(ie) = iet_('eQTQT');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E29';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTRC';
        ielftype(ie) = iet_('eQTRC');
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E30';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTRC';
        ielftype(ie) = iet_('eQTRC');
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E31';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTQT';
        ielftype(ie) = iet_('eQTQT');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E32';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTQT';
        ielftype(ie) = iet_('eQTQT');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E33';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E34';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eINV';
        ielftype(ie) = iet_('eINV');
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E35';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E36';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQTRC';
        ielftype(ie) = iet_('eQTRC');
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E37';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eINV';
        ielftype(ie) = iet_('eINV');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E38';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PRRC';
        ielftype(ie) = iet_('en2PRRC');
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E39';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E40';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E41';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E42';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E43';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E44';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E45';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQT';
        ielftype(ie) = iet_('eQT');
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,0.9,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = -1.231060;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = -1.231060;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        ig = ig_('C0');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = -1.231060;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = -1.231060;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.231060;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.034750;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        pbm.grelw{ig}(posel) = -0.00975;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.034750;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        pbm.grelw{ig}(posel) = -0.00975;
        ig = ig_('C3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.034750;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E11');
        pbm.grelw{ig}(posel) = -0.00975;
        ig = ig_('C4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E12');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.034750;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E13');
        pbm.grelw{ig}(posel) = -0.00975;
        ig = ig_('C5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.034750;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E15');
        pbm.grelw{ig}(posel) = -0.00975;
        ig = ig_('C6');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E16');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E17');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C7');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E19');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E20');
        pbm.grelw{ig}(posel) = 0.002;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E21');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.002;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E22');
        pbm.grelw{ig}(posel) = -0.002;
        ig = ig_('C8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E23');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.002;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E24');
        pbm.grelw{ig}(posel) = 0.002;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.002;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E26');
        pbm.grelw{ig}(posel) = -0.002;
        ig = ig_('C9');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E27');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E28');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E29');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 500.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E30');
        pbm.grelw{ig}(posel) = -500.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E31');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_('C10');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E32');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E33');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E34');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 500.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E35');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E36');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -500.0;
        ig = ig_('C11');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E37');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.9;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E38');
        pbm.grelw{ig}(posel) = -0.002;
        ig = ig_('C13');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E39');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C14');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E40');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C15');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E41');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C16');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E42');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C17');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E43');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C18');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E44');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('C19');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E45');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               174.788807
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MN-16-20';
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

    case 'eINV'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 1.0/EV_(1);
        if(nargout>1)
            g_(1,1) = -1.0/EV_(1)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0/EV_(1)^3;
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

    case 'eQT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)/EV_(2);
        if(nargout>1)
            g_(1,1) = 1/EV_(2);
            g_(2,1) = -EV_(1)/EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -1.0/EV_(2)^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = (2.0*EV_(1))/EV_(2)^3;
                varargout{3} = H_;
            end
        end

    case 'eQTQT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XW = EV_(1)*EV_(3);
        YZ = EV_(2)*EV_(4);
        varargout{1} = XW/YZ;
        if(nargout>1)
            g_(1,1) = EV_(3)/YZ;
            g_(2,1) = -XW/(EV_(2)*YZ);
            g_(3,1) = EV_(1)/YZ;
            g_(4,1) = -XW/(YZ*EV_(4));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = -EV_(3)/(EV_(2)*YZ);
                H_(2,1) = H_(1,2);
                H_(1,3) = 1.0/YZ;
                H_(3,1) = H_(1,3);
                H_(1,4) = -EV_(3)/(EV_(4)*YZ);
                H_(4,1) = H_(1,4);
                H_(2,2) = (2.0*XW)/(EV_(2)^2*YZ);
                H_(2,3) = -EV_(1)/(EV_(2)*YZ);
                H_(3,2) = H_(2,3);
                H_(2,4) = XW/YZ^2;
                H_(4,2) = H_(2,4);
                H_(3,4) = -EV_(1)/(EV_(4)*YZ);
                H_(4,3) = H_(3,4);
                H_(4,4) = (2.0*XW)/(EV_(4)^2*YZ);
                varargout{3} = H_;
            end
        end

    case 'en2PRRC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)/EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)/EV_(3);
            g_(2,1) = EV_(1)/EV_(3);
            g_(3,1) = -EV_(1)*EV_(2)/EV_(3)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = 1.0/EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = -EV_(2)/EV_(3)^2;
                H_(3,1) = H_(1,3);
                H_(2,3) = -EV_(1)/EV_(3)^2;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*EV_(1)*EV_(2)/EV_(3)^3;
                varargout{3} = H_;
            end
        end

    case 'eQTRC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        YZ = EV_(2)*EV_(3);
        varargout{1} = EV_(1)/YZ;
        if(nargout>1)
            g_(1,1) = 1.0/YZ;
            g_(2,1) = -EV_(1)/(EV_(2)*YZ);
            g_(3,1) = -EV_(1)/(EV_(3)*YZ);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = -1.0/(EV_(2)*YZ);
                H_(2,1) = H_(1,2);
                H_(1,3) = -1.0/(EV_(3)*YZ);
                H_(3,1) = H_(1,3);
                H_(2,2) = (2.0*EV_(1))/(EV_(2)^2*YZ);
                H_(2,3) = EV_(1)/YZ^2;
                H_(3,2) = H_(2,3);
                H_(3,3) = (2.0*EV_(1))/(EV_(3)^2*YZ);
                varargout{3} = H_;
            end
        end

    case 'eSQQT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^2/EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)/EV_(2);
            g_(2,1) = -EV_(1)^2/EV_(2)^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0/EV_(2);
                H_(1,2) = -2.0*EV_(1)/EV_(2)^2;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1)^2/EV_(2)^3;
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

