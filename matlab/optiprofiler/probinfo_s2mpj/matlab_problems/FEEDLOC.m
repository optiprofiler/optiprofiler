function varargout = FEEDLOC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : FEEDLOC
%    *********
% 
%    Feed tray location & determination of optimum number of trays 
%    in a distillation column
% 
%    SIF input: S. Leyffer, October 1997
% 
%    classification = 'C-CLOR2-AN-90-259'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'FEEDLOC';

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
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('M') = 2;
        v_('NMAX') = 12;
        v_('NMAX-1') = -1+v_('NMAX');
        v_('F') = 100.0;
        v_('AL1') = 1.0;
        v_('AL2') = 5.13435;
        v_('XF1') = 0.80;
        v_('XF2') = 0.20;
        v_('SPEC') = 0.001;
        v_('BIGM') = 1000.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('NMAX')
            [iv,ix_] = s2mpjlib('ii',['S',int2str(I)],ix_);
            pb.xnames{iv} = ['S',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['W',int2str(I)],ix_);
            pb.xnames{iv} = ['W',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','N',ix_);
        pb.xnames{iv} = 'N';
        for I=v_('1'):v_('NMAX')
            for J=v_('1'):v_('M')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Y',int2str(I),',',int2str(J)];
            end
        end
        [iv,ix_] = s2mpjlib('ii','L',ix_);
        pb.xnames{iv} = 'L';
        [iv,ix_] = s2mpjlib('ii','V',ix_);
        pb.xnames{iv} = 'V';
        [iv,ix_] = s2mpjlib('ii','R',ix_);
        pb.xnames{iv} = 'R';
        [iv,ix_] = s2mpjlib('ii','P1',ix_);
        pb.xnames{iv} = 'P1';
        [iv,ix_] = s2mpjlib('ii','P2',ix_);
        pb.xnames{iv} = 'P2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('R');
        valA(end+1) = 1.0;
        for I=v_('1'):v_('NMAX')
            [ig,ig_] = s2mpjlib('ii','FENTR',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'FENTR';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii','NTRAY',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'NTRAY';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii','NDEF1',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'NDEF1';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = 1.0;
            v_('RI') = I;
            [ig,ig_] = s2mpjlib('ii','NDEF2',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'NDEF2';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = v_('RI');
        end
        [ig,ig_] = s2mpjlib('ii','NDEF1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NDEF1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('N');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','NDEF2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NDEF2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('N');
        valA(end+1) = -1.0;
        for I=v_('1'):v_('NMAX-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['NIL',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['NIL',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(round(v_('I+1')))]);
            valA(end+1) = -1.0;
        end
        for I=v_('1'):v_('NMAX')
            v_('RI') = I;
            [ig,ig_] = s2mpjlib('ii','ENTX',ig_);
            gtype{ig}  = '<=';
            cnames{ig} = 'ENTX';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = v_('RI');
            v_('RI') = -1.0*v_('RI');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = v_('RI');
        end
        for I=v_('1'):v_('NMAX')
            [ig,ig_] = s2mpjlib('ii',['LASTX',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['LASTX',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for I=v_('1'):v_('NMAX')
            [ig,ig_] = s2mpjlib('ii',['ZNOT',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['ZNOT',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = 1.0;
            for K=I:v_('NMAX')
                [ig,ig_] = s2mpjlib('ii',['ZNOT',int2str(I)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['ZNOT',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['S',int2str(K)]);
                valA(end+1) = -1.0;
            end
        end
        for I=v_('1'):v_('NMAX')
            [ig,ig_] = s2mpjlib('ii',['FEEDX',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['FEEDX',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for I=v_('2'):v_('NMAX')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['WNES1u',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['WNES1u',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['WNES2u',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['WNES2u',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(round(v_('I-1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = 1.0;
        end
        for I=v_('1'):v_('NMAX')
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['PE1',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['PE1',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Y',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['PE2',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['PE2',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Y',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['PE3',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['PE3',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                [ig,ig_] = s2mpjlib('ii',['PE4',int2str(I)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['PE4',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['PE1',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['PE1',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['PE2',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['PE2',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['PE3',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['PE3',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['PE4',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['PE4',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = -1.0;
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['XNOT',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['XNOT',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Z',int2str(I)]);
                valA(end+1) = -1.0;
                [ig,ig_] = s2mpjlib('ii',['YNOT',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['YNOT',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Y',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['Z',int2str(I)]);
                valA(end+1) = -1.0;
            end
        end
        for I=v_('1'):v_('NMAX')
            v_('TEMP') = -1.0*v_('AL1');
            [ig,ig_] = s2mpjlib('ii',['PHEE',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PHEE',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('1')))]);
            valA(end+1) = v_('TEMP');
        end
        [ig,ig_] = s2mpjlib('ii','DEFL',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'DEFL';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('L');
        valA(end+1) = 1.0;
        for J=v_('1'):v_('M')
            v_('TEMP') = -1.0*v_('F');
            [ig,ig_] = s2mpjlib('ii',['CMB1u',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CMB1u',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('2'))),',',int2str(J)]);
            valA(end+1) = v_('TEMP');
        end
        for I=v_('2'):v_('NMAX')
            for J=v_('1'):v_('M')
                v_('TEMP') = -1.0*v_('BIGM');
                [ig,ig_] = s2mpjlib('ii',['CMBN1',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['CMBN1',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['S',int2str(I)]);
                valA(end+1) = v_('BIGM');
                [ig,ig_] = s2mpjlib('ii',['CMBN2',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['CMBN2',int2str(I),',',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['S',int2str(I)]);
                valA(end+1) = v_('TEMP');
            end
        end
        for I=v_('2'):v_('NMAX-1')
            v_('TEMP1') = v_('F')*v_(['XF',int2str(round(v_('M')))]);
            v_('TEMP1') = -1.0*v_('TEMP1');
            [ig,ig_] = s2mpjlib('ii',['CMB1',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['CMB1',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = v_('TEMP');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = v_('BIGM');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = v_('TEMP1');
            [ig,ig_] = s2mpjlib('ii',['CMB2',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CMB2',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = v_('BIGM');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['Z',int2str(I)]);
            valA(end+1) = v_('TEMP');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['W',int2str(I)]);
            valA(end+1) = v_('TEMP1');
        end
        for I=v_('3'):v_('NMAX')
            [ig,ig_] = s2mpjlib('ii',['RECR',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['RECR',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['S',int2str(I)]);
            valA(end+1) = v_('BIGM');
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
        pbm.gconst(ig_('FENTR')) = 1.0;
        pbm.gconst(ig_('NTRAY')) = 1.0;
        for I=v_('2'):v_('NMAX')
            pbm.gconst(ig_(['WNES1u',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['WNES2u',int2str(I)])) = 1.0;
        end
        v_('TEMP') = -1.0*v_('BIGM');
        for I=v_('2'):v_('NMAX')
            for J=v_('1'):v_('M')
                pbm.gconst(ig_(['CMBN1',int2str(I),',',int2str(J)])) = v_('BIGM');
                pbm.gconst(ig_(['CMBN2',int2str(I),',',int2str(J)])) = v_('TEMP');
            end
        end
        for I=v_('2'):v_('NMAX-1')
            pbm.gconst(ig_(['CMB1',int2str(I)])) = v_('BIGM');
            v_('TEMP') = -1.0*v_('BIGM');
            pbm.gconst(ig_(['CMB2',int2str(I)])) = v_('TEMP');
        end
        v_('TEMP') = v_(['XF',int2str(round(v_('1')))])*v_('SPEC');
        v_('TEMP1') = v_('TEMP')*v_('F');
        v_('RHS') = v_('TEMP1')+v_('BIGM');
        for I=v_('3'):v_('NMAX')
            pbm.gconst(ig_(['RECR',int2str(I)])) = v_('RHS');
        end
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        for I=v_('1'):v_('NMAX')
            grange(ig_(['PE1',int2str(I)])) = 2.0;
            grange(ig_(['PE2',int2str(I)])) = 2.0;
            grange(ig_(['PE3',int2str(I)])) = 2.0;
            grange(ig_(['PE4',int2str(I)])) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('NMAX')
            pb.xupper(ix_(['Z',int2str(I)])) = 1.0;
            pb.xupper(ix_(['W',int2str(I)])) = 1.0;
            pb.xupper(ix_(['S',int2str(I)])) = 1.0;
            for J=v_('1'):v_('M')
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)])) = 1.0;
                pb.xupper(ix_(['Y',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        pb.xlower(ix_('N'),1) = 3.0;
        v_('TEMP') = v_('NMAX');
        pb.xupper(ix_('N')) = v_('TEMP');
        pb.xlower(ix_('P2'),1) = 80.0;
        pb.xupper(ix_('P2'),1) = 80.0;
        pb.xupper(ix_('L')) = v_('F');
        pb.xupper(ix_('V')) = v_('F');
        pb.xupper(ix_('P1')) = v_('F');
        pb.xupper(ix_('R')) = 5.0;
        pb.xlower(ix_(['W',int2str(round(v_('1')))]),1) = 0.0;
        pb.xupper(ix_(['W',int2str(round(v_('1')))]),1) = 0.0;
        pb.xlower(ix_(['W',int2str(round(v_('2')))]),1) = 0.0;
        pb.xupper(ix_(['W',int2str(round(v_('2')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        pb.y0 = 0.5*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eA2PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('NMAX')
            for K=v_('1'):v_('M')
                ename = ['PHE',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = ['X',int2str(I),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['AL',int2str(K)]);
            end
        end
        ename = 'DEFLE';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA2PROD';
        ielftype(ie) = iet_('eA2PROD');
        vname = 'R';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'P1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('A',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        for J=v_('1'):v_('M')
            ename = ['CMB11u',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'P2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            ename = ['CMB12u',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'V';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(round(v_('1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            ename = ['CMB13u',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('2'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
        end
        for I=v_('2'):v_('NMAX')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('M')
                ename = ['CM11',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'L';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                ename = ['CM12',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'P1';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                ename = ['CM13',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'V';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
                ename = ['CM21',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'L';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                ename = ['CM22',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'P1';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = 1.0;
                ename = ['CM23',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = 'V';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = -1.0;
            end
        end
        for I=v_('2'):v_('NMAX-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            ename = ['C11',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            for K=I:v_('NMAX')
                ename = ['C12',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = ['X',int2str(I),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['W',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('F');
            end
            ename = ['C13',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'V';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            ename = ['C14',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
            for K=v_('I+1'):v_('NMAX')
                v_('TEMP') = -1.0*v_('F');
                ename = ['C15',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['W',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('TEMP');
            end
            ename = ['C16',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'V';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
            ename = ['C21',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            for K=v_('1'):v_('NMAX')
                ename = ['C22',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = ['X',int2str(I),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['W',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('F');
            end
            ename = ['C23',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'V';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
            ename = ['C24',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
            for K=v_('I+1'):v_('NMAX')
                v_('TEMP') = -1.0*v_('F');
                ename = ['C25',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eA2PROD';
                ielftype(ie) = iet_('eA2PROD');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['W',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('A',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('TEMP');
            end
            ename = ['C26',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'V';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
        end
        for I=v_('3'):v_('NMAX')
            ename = ['REC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eA2PROD';
            ielftype(ie) = iet_('eA2PROD');
            vname = 'P1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I),',',int2str(round(v_('1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.5);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NMAX')
            for K=v_('1'):v_('M')
                ig = ig_(['PHEE',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['PHE',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        ig = ig_('DEFL');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('DEFLE');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        for J=v_('1'):v_('M')
            ig = ig_(['CMB1u',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CMB11u',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['CMB12u',int2str(J)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CMB13u',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for I=v_('2'):v_('NMAX')
            for J=v_('1'):v_('M')
                ig = ig_(['CMBN1',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['CM11',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['CM12',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['CM13',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['CMBN2',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['CM21',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['CM22',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['CM23',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for I=v_('2'):v_('NMAX-1')
            ig = ig_(['CMB1',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C11',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C13',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C14',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C16',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['CMB2',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C21',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C23',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C24',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['C26',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            for K=I:v_('NMAX')
                ig = ig_(['CMB1',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C12',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['CMB2',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C22',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
            v_('I+1') = 1+I;
            for K=v_('I+1'):v_('NMAX')
                ig = ig_(['CMB1',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C15',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['CMB2',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['C25',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for I=v_('3'):v_('NMAX')
            ig = ig_(['RECR',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['REC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = grange(gegrps);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-AN-90-259';
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

    case 'eA2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(2);
            g_(2,1) = pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
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

