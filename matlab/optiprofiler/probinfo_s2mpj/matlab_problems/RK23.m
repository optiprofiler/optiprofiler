function varargout = RK23(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : RK23
%    *********
% 
%    Find coefficients for an embedded pair of explicit 2nd
%    and 3rd order Runge Kutta Method such that the leading
%    term in the local truncation error is minimized.
% 
%    Source:
%    Similar ideas for 4th and 5th order pairs are discussed in:
%    Hairer, Norsett and Wanner, Solving Ordinary Differential
%    Equations I, Springer 1980, page 158 ff.
% 
%    SIF input: S. Leyffer, January 1997.
% 
%    classification = 'C-CLOR2-RN-17-11'
% 
% 
%    ... COMPUTED PARAMETERS
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'RK23';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('ONE') = 1.0;
        v_('THREE') = 3.0;
        v_('FOUR') = 4.0;
        v_('SIX') = 6.0;
        v_('ONETHIRD') = v_('ONE')/v_('THREE');
        v_('ONESIXTH') = v_('ONE')/v_('SIX');
        v_('FOURSIXTH') = v_('FOUR')/v_('SIX');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','C2',ix_);
        pb.xnames{iv} = 'C2';
        [iv,ix_] = s2mpjlib('ii','A21',ix_);
        pb.xnames{iv} = 'A21';
        [iv,ix_] = s2mpjlib('ii','C3',ix_);
        pb.xnames{iv} = 'C3';
        [iv,ix_] = s2mpjlib('ii','A31',ix_);
        pb.xnames{iv} = 'A31';
        [iv,ix_] = s2mpjlib('ii','A32',ix_);
        pb.xnames{iv} = 'A32';
        [iv,ix_] = s2mpjlib('ii','B1',ix_);
        pb.xnames{iv} = 'B1';
        [iv,ix_] = s2mpjlib('ii','B2',ix_);
        pb.xnames{iv} = 'B2';
        [iv,ix_] = s2mpjlib('ii','B3',ix_);
        pb.xnames{iv} = 'B3';
        [iv,ix_] = s2mpjlib('ii','BB1',ix_);
        pb.xnames{iv} = 'BB1';
        [iv,ix_] = s2mpjlib('ii','BB2',ix_);
        pb.xnames{iv} = 'BB2';
        [iv,ix_] = s2mpjlib('ii','BB3',ix_);
        pb.xnames{iv} = 'BB3';
        [iv,ix_] = s2mpjlib('ii','TP1',ix_);
        pb.xnames{iv} = 'TP1';
        [iv,ix_] = s2mpjlib('ii','TM1',ix_);
        pb.xnames{iv} = 'TM1';
        [iv,ix_] = s2mpjlib('ii','TP2',ix_);
        pb.xnames{iv} = 'TP2';
        [iv,ix_] = s2mpjlib('ii','TM2',ix_);
        pb.xnames{iv} = 'TM2';
        [iv,ix_] = s2mpjlib('ii','TP3',ix_);
        pb.xnames{iv} = 'TP3';
        [iv,ix_] = s2mpjlib('ii','TM3',ix_);
        pb.xnames{iv} = 'TM3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('TP1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TP2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','ROWS1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ROWS1';
        iv = ix_('A21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('C2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','ROWS2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ROWS2';
        iv = ix_('A31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('A32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('C3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FIRST2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'FIRST2';
        iv = ix_('B1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('B2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('B3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','FIRST3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'FIRST3';
        iv = ix_('BB1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('BB2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('BB3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','SECND2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SECND2';
        [ig,ig_] = s2mpjlib('ii','SECND3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SECND3';
        [ig,ig_] = s2mpjlib('ii','THIRD31',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'THIRD31';
        [ig,ig_] = s2mpjlib('ii','THIRD32',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'THIRD32';
        [ig,ig_] = s2mpjlib('ii','ART1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ART1';
        iv = ix_('TP1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','ART2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ART2';
        iv = ix_('TP2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','ART3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ART3';
        iv = ix_('TP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('TM3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
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
        pbm.gconst(ig_('ROWS1')) = 0.0;
        pbm.gconst(ig_('ROWS2')) = 0.0;
        pbm.gconst(ig_('FIRST2')) = 1.0;
        pbm.gconst(ig_('FIRST3')) = 1.0;
        pbm.gconst(ig_('SECND2')) = 0.5;
        pbm.gconst(ig_('SECND3')) = 0.5;
        pbm.gconst(ig_('THIRD31')) = v_('ONETHIRD');
        pbm.gconst(ig_('THIRD32')) = v_('ONESIXTH');
        pbm.gconst(ig_('ART1')) = 1.0;
        pbm.gconst(ig_('ART2')) = 1.0;
        pbm.gconst(ig_('ART3')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower(ix_('TP1'),1) = 0.0;
        pb.xlower(ix_('TM1'),1) = 0.0;
        pb.xlower(ix_('TP2'),1) = 0.0;
        pb.xlower(ix_('TM2'),1) = 0.0;
        pb.xlower(ix_('TP3'),1) = 0.0;
        pb.xlower(ix_('TM3'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'C2'))
            pb.x0(ix_('C2'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('C2')),1) = 1.0;
        end
        if(isKey(ix_,'A21'))
            pb.x0(ix_('A21'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('A21')),1) = 1.0;
        end
        if(isKey(ix_,'C3'))
            pb.x0(ix_('C3'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('C3')),1) = 0.5;
        end
        if(isKey(ix_,'A31'))
            pb.x0(ix_('A31'),1) = 0.25;
        else
            pb.y0(find(pbm.congrps==ig_('A31')),1) = 0.25;
        end
        if(isKey(ix_,'A32'))
            pb.x0(ix_('A32'),1) = 0.25;
        else
            pb.y0(find(pbm.congrps==ig_('A32')),1) = 0.25;
        end
        if(isKey(ix_,'B1'))
            pb.x0(ix_('B1'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 0.5;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = 0.5;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = 0.5;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = 0.0;
        end
        pb.x0(ix_('BB1'),1) = v_('ONESIXTH');
        pb.x0(ix_('BB2'),1) = v_('ONESIXTH');
        pb.x0(ix_('BB3'),1) = v_('FOURSIXTH');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'ePRODS',iet_);
        elftv{it}{1} = 'W1';
        elftv{it}{2} = 'W2';
        [it,iet_] = s2mpjlib( 'ii', 'ePRODQ',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        [it,iet_] = s2mpjlib( 'ii', 'eTPROD',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        [it,iet_] = s2mpjlib( 'ii', 'eQPROD',iet_);
        elftv{it}{1} = 'Z1';
        elftv{it}{2} = 'Z2';
        elftv{it}{3} = 'Z3';
        elftv{it}{4} = 'Z4';
        [it,iet_] = s2mpjlib( 'ii', 'eTPRODS',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'B2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'B3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'BB2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePROD';
        ielftype(ie) = iet_('ePROD');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODS';
        ielftype(ie) = iet_('ePRODS');
        vname = 'BB2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODS';
        ielftype(ie) = iet_('ePRODS');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eTPROD';
        ielftype(ie) = iet_('eTPROD');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'A32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODQ';
        ielftype(ie) = iet_('ePRODQ');
        vname = 'BB2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'ePRODQ';
        ielftype(ie) = iet_('ePRODQ');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eQPROD';
        ielftype(ie) = iet_('eQPROD');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'A32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eTPRODS';
        ielftype(ie) = iet_('eTPRODS');
        vname = 'BB3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'A32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'C2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('SECND2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('SECND3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('THIRD31');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('THIRD32');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ART1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 4.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        pbm.grelw{ig}(posel) = 4.0;
        ig = ig_('ART2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 8.0;
        ig = ig_('ART3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 12.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-RN-17-11';
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
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'ePRODS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*(EV_(2)^2.0);
        if(nargout>1)
            g_(1,1) = EV_(2)^2;
            g_(2,1) = EV_(1)*2.0*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 2.0*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1);
                varargout{3} = H_;
            end
        end

    case 'ePRODQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*(EV_(2)^3.0);
        if(nargout>1)
            g_(1,1) = EV_(2)^3.0;
            g_(2,1) = EV_(1)*3.0*(EV_(2)^2.0);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 3.0*(EV_(2)^2.0);
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)*6.0*EV_(2);
                varargout{3} = H_;
            end
        end

    case 'eTPROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    case 'eQPROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3)*EV_(4);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3)*EV_(4);
            g_(2,1) = EV_(1)*EV_(3)*EV_(4);
            g_(3,1) = EV_(1)*EV_(2)*EV_(4);
            g_(4,1) = EV_(1)*EV_(2)*EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = EV_(3)*EV_(4);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*EV_(4);
                H_(3,1) = H_(1,3);
                H_(1,4) = EV_(2)*EV_(3);
                H_(4,1) = H_(1,4);
                H_(2,3) = EV_(1)*EV_(4);
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*EV_(3);
                H_(4,2) = H_(2,4);
                H_(3,4) = EV_(1)*EV_(2);
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    case 'eTPRODS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*(EV_(3)^2.0);
        if(nargout>1)
            g_(1,1) = EV_(2)*(EV_(3)^2.0);
            g_(2,1) = EV_(1)*(EV_(3)^2.0);
            g_(3,1) = EV_(1)*EV_(2)*2.0*EV_(3);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3)^2.0;
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2)*2.0*EV_(3);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1)*2.0*EV_(3);
                H_(3,2) = H_(2,3);
                H_(3,3) = EV_(1)*EV_(2)*2.0;
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

