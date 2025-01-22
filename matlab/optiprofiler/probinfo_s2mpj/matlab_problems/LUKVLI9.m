function varargout = LUKVLI9(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKVLI9
%    *********
% 
%    Source: Problem 5.9, the modified Brown function with
%    simplified seven-diagonal constraints, due to L. Luksan and J. Vlcek,
%    "Sparse and partially separable test problems for 
%    unconstrained and equality constrained optimization",
%    Technical Report 767, Inst. Computer Science, Academy of Sciences
%    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
% 
%    Equality constraints changed to inequalities
% 
%    SIF input: Nick Gould, April 2001
% 
%    classification = 'C-COOR2-AY-V-V'
% 
%    some useful parameters, including N, the number of variables.
% 
%       Alternative values for the SIF file parameters:
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   10000          $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKVLI9';

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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   100000         $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('N-3') = -3+v_('N');
        v_('N-4') = -4+v_('N');
        v_('N-5') = -5+v_('N');
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
        for I=v_('1'):v_('N/2')
            v_('2I') = 2*I;
            v_('2I-1') = -1+v_('2I');
            [ig,ig_] = s2mpjlib('ii',['OBJ1',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I-1')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii','OBJ2',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I-1')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['OBJ3',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I-1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2I')))]);
            valA(end+1) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('1')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = 4.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('2')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('1')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('3')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('2')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('2')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('2')))]);
        valA(end+1) = 6.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('3')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('2')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('2')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('4')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('3')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('3')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('3')))]);
        valA(end+1) = 6.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('4')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('3')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('3')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('5')))]);
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('4')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-2')))]);
        valA(end+1) = 6.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-1')))]);
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('4')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-4')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('4')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('4')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-5')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('5')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('5')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-1')))]);
        valA(end+1) = 6.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-3')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('5')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('5')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-4')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('6')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('6')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
        valA(end+1) = 2.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-3')))]);
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('6')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['C',int2str(round(v_('6')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['X',int2str(round(v_('N-2')))]);
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
        pbm.gconst(ig_(['C',int2str(round(v_('2')))])) = 2.0;
        pbm.gconst(ig_(['C',int2str(round(v_('3')))])) = 2.0;
        pbm.gconst(ig_(['C',int2str(round(v_('4')))])) = 2.0;
        pbm.gconst(ig_(['C',int2str(round(v_('5')))])) = 2.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCUBEP',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = ['C1',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C1',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C4',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C4',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C5',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C5',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('5')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C4',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C4',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('5')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C5',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C5',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C6',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C6',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('6')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C4',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C4',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C5',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C5',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-4')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C4',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C4',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-3')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCUBEP';
        ielftype(ie) = iet_('eCUBEP');
        ename = ['C1',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C1',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C2',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C2',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-1')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['C3',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        ename = ['C3',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = ['X',int2str(round(v_('N-2')))];
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gAL2',igt_);
        [it,igt_] = s2mpjlib('ii','gEXP20',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N/2')
            ig = ig_(['OBJ1',int2str(I)]);
            pbm.grftype{ig} = 'gAL2';
            ig = ig_(['OBJ3',int2str(I)]);
            pbm.grftype{ig} = 'gEXP20';
        end
        ig = ig_(['C',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C4',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C5',int2str(round(v_('2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C4',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C5',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C6',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        ig = ig_(['C',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C4',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C5',int2str(round(v_('4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_(['C',int2str(round(v_('5')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -4.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C4',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_(['C',int2str(round(v_('6')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C1',int2str(round(v_('6')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C2',int2str(round(v_('6')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['C3',int2str(round(v_('6')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               9.98933E+01
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AY-V-V';
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

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eCUBEP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^3-EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 3.0*EV_(1)^2-EV_(2);
            g_(2,1) = -EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 6.0*EV_(1);
                H_(1,2) = -1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gAL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 0.001*GVAR_*GVAR_;
        if(nargout>1)
            g_ = 0.002*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 0.002;
                varargout{3} = H_;
            end
        end

    case 'gEXP20'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        EXP20A = exp(20.0*GVAR_);
        varargout{1} = EXP20A;
        if(nargout>1)
            g_ = 20.0*EXP20A;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 400.0*EXP20A;
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

