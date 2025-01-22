function varargout = LOADBAL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    The problem arises in the field of computer networks and parallel
%    computation.  It deals with the static load balancing in a tree
%    computer network with two-way traffic.  A set of heterogeneous host
%    computers are interconnected, in which each node processes jobs (the 
%    jobs arriving at each node according to a time invariant Poisson process) 
%    locally or sends it to a remote node,.  In the latter case, there is a
%    communication delay of forwarding the job and getting a response back.
%    The problem is then to minimize the mean response time of a job.
% 
%    The example considered here features 11 computers arranged as follows:
% 
%          1      6      9
%           \     |     /
%            \    |    /
%         2---4---5---8---10
%            /    |    \
%           /     |     \
%          3      7      11
% 
%    Source:
%    J. Li and H. Kameda,
%    "Optimal load balancing in tree network with two-way traffic",
%    Computer networks and ISDN systems, vol. 25, pp. 1335-1348, 1993.
% 
%    SIF input: Masha Sosonkina, Virginia Tech., 1995.
% 
%    classification = 'C-COLR2-MN-31-31'
% 
%  Parameter assignment.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LOADBAL';

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
        v_('P1') = 3;
        v_('N') = 11;
        v_('NLINK') = 20;
        v_('NLINK-3') = 17;
        v_('NLINK-4') = 16;
        v_('4C') = 4;
        v_('5C') = 5;
        v_('6C') = 6;
        v_('7C') = 7;
        v_('8C') = 8;
        v_('FI') = 514.0;
        v_('0.2*FI') = 0.2*v_('FI');
        v_('0.0125*FI') = 0.0125*v_('FI');
        v_('0.05*FI') = 0.05*v_('FI');
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('0.2*FI');
        end
        for I=v_('1'):v_('NLINK')
            [ig,ig_] = s2mpjlib('ii',['CNST',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['CNST',int2str(I)];
            [ig,ig_] = s2mpjlib('ii',['GA',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('0.0125*FI');
            [ig,ig_] = s2mpjlib('ii',['GB',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('0.05*FI');
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['N',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['N',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = ig_.Count;
        [iv,ix_] = s2mpjlib('ii','X4,1',ix_);
        pb.xnames{iv} = 'X4,1';
        icA(end+1) = iv;
        irA(end+1) = ig_('N1');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X1,4',ix_);
        pb.xnames{iv} = 'X1,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('N1');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X4,2',ix_);
        pb.xnames{iv} = 'X4,2';
        icA(end+1) = iv;
        irA(end+1) = ig_('N2');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X2,4',ix_);
        pb.xnames{iv} = 'X2,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('N2');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X4,3',ix_);
        pb.xnames{iv} = 'X4,3';
        icA(end+1) = iv;
        irA(end+1) = ig_('N3');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X3,4',ix_);
        pb.xnames{iv} = 'X3,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('N3');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X4,5',ix_);
        pb.xnames{iv} = 'X4,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X5,4',ix_);
        pb.xnames{iv} = 'X5,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N4');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X5,6',ix_);
        pb.xnames{iv} = 'X5,6';
        icA(end+1) = iv;
        irA(end+1) = ig_('N6');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X6,5',ix_);
        pb.xnames{iv} = 'X6,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('N6');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X5,7',ix_);
        pb.xnames{iv} = 'X5,7';
        icA(end+1) = iv;
        irA(end+1) = ig_('N7');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X7,5',ix_);
        pb.xnames{iv} = 'X7,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('N7');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X5,8',ix_);
        pb.xnames{iv} = 'X5,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X8,5',ix_);
        pb.xnames{iv} = 'X8,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N5');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X8,9',ix_);
        pb.xnames{iv} = 'X8,9';
        icA(end+1) = iv;
        irA(end+1) = ig_('N9');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X9,8',ix_);
        pb.xnames{iv} = 'X9,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('N9');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X8,10',ix_);
        pb.xnames{iv} = 'X8,10';
        icA(end+1) = iv;
        irA(end+1) = ig_('N10');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X10,8',ix_);
        pb.xnames{iv} = 'X10,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('N10');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X8,11',ix_);
        pb.xnames{iv} = 'X8,11';
        icA(end+1) = iv;
        irA(end+1) = ig_('N11');
        valA(end+1) = 1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = -1.00;
        [iv,ix_] = s2mpjlib('ii','X11,8',ix_);
        pb.xnames{iv} = 'X11,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('N11');
        valA(end+1) = -1.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('N8');
        valA(end+1) = 1.00;
        [iv,ix_] = s2mpjlib('ii','X4,1',ix_);
        pb.xnames{iv} = 'X4,1';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST1');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST2');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X1,4',ix_);
        pb.xnames{iv} = 'X1,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST1');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST2');
        valA(end+1) = 20.00;
        [iv,ix_] = s2mpjlib('ii','X4,2',ix_);
        pb.xnames{iv} = 'X4,2';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST3');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST4');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X2,4',ix_);
        pb.xnames{iv} = 'X2,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST3');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST4');
        valA(end+1) = 20.0;
        [iv,ix_] = s2mpjlib('ii','X4,3',ix_);
        pb.xnames{iv} = 'X4,3';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST5');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST6');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X3,4',ix_);
        pb.xnames{iv} = 'X3,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST5');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST6');
        valA(end+1) = 20.00;
        [iv,ix_] = s2mpjlib('ii','X5,6',ix_);
        pb.xnames{iv} = 'X5,6';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST7');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST8');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X6,5',ix_);
        pb.xnames{iv} = 'X6,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST7');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST8');
        valA(end+1) = 20.0;
        [iv,ix_] = s2mpjlib('ii','X5,7',ix_);
        pb.xnames{iv} = 'X5,7';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST9');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST10');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X7,5',ix_);
        pb.xnames{iv} = 'X7,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST9');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST10');
        valA(end+1) = 20.0;
        [iv,ix_] = s2mpjlib('ii','X8,9',ix_);
        pb.xnames{iv} = 'X8,9';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST11');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST12');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X9,8',ix_);
        pb.xnames{iv} = 'X9,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST11');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST12');
        valA(end+1) = 20.0;
        [iv,ix_] = s2mpjlib('ii','X8,10',ix_);
        pb.xnames{iv} = 'X8,10';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST13');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST14');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X10,8',ix_);
        pb.xnames{iv} = 'X10,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST13');
        valA(end+1) = 80.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST14');
        valA(end+1) = 20.0;
        [iv,ix_] = s2mpjlib('ii','X8,11',ix_);
        pb.xnames{iv} = 'X8,11';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST15');
        valA(end+1) = 20.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST16');
        valA(end+1) = 80.00;
        [iv,ix_] = s2mpjlib('ii','X11,8',ix_);
        pb.xnames{iv} = 'X11,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST15');
        valA(end+1) = 80.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST16');
        valA(end+1) = 20.00;
        [iv,ix_] = s2mpjlib('ii','X4,5',ix_);
        pb.xnames{iv} = 'X4,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST17');
        valA(end+1) = 20.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST18');
        valA(end+1) = 80.0;
        [iv,ix_] = s2mpjlib('ii','X5,4',ix_);
        pb.xnames{iv} = 'X5,4';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST17');
        valA(end+1) = 80.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST18');
        valA(end+1) = 20.00;
        [iv,ix_] = s2mpjlib('ii','X5,8',ix_);
        pb.xnames{iv} = 'X5,8';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST19');
        valA(end+1) = 20.00;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST20');
        valA(end+1) = 80.0;
        [iv,ix_] = s2mpjlib('ii','X8,5',ix_);
        pb.xnames{iv} = 'X8,5';
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST19');
        valA(end+1) = 80.0;
        icA(end+1) = iv;
        irA(end+1) = ig_('CNST20');
        valA(end+1) = 20.00;
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
            icA(end+1) = iv;
            irA(end+1) = ig_(['N',int2str(I)]);
            valA(end+1) = -1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
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
        pbm.gconst(ig_('N1')) = -95.0;
        pbm.gconst(ig_('N2')) = -95.0;
        pbm.gconst(ig_('N3')) = -19.0;
        pbm.gconst(ig_('N4')) = -70.0;
        pbm.gconst(ig_('N5')) = -70.0;
        pbm.gconst(ig_('N6')) = -19.0;
        pbm.gconst(ig_('N7')) = -19.0;
        pbm.gconst(ig_('N8')) = -70.0;
        pbm.gconst(ig_('N9')) = -19.0;
        pbm.gconst(ig_('N10')) = -19.0;
        pbm.gconst(ig_('N11')) = -19.0;
        v_('CIJE') = 999.99;
        for I=v_('1'):v_('NLINK-4')
            pbm.gconst(ig_(['CNST',int2str(I)])) = v_('CIJE');
        end
        v_('CIJE') = 9999.99;
        for I=v_('NLINK-3'):v_('NLINK')
            pbm.gconst(ig_(['CNST',int2str(I)])) = v_('CIJE');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xupper(ix_('B1')) = 99.99;
        pb.xupper(ix_('B2')) = 99.99;
        pb.xupper(ix_('B4')) = 99.99;
        pb.xupper(ix_('B5')) = 99.99;
        pb.xupper(ix_('B8')) = 99.99;
        pb.xupper(ix_('B3')) = 19.99;
        pb.xupper(ix_('B6')) = 19.99;
        pb.xupper(ix_('B7')) = 19.99;
        pb.xupper(ix_('B9')) = 19.99;
        pb.xupper(ix_('B10')) = 19.99;
        pb.xupper(ix_('B11')) = 19.99;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1,4'))
            pb.x0(ix_('X1,4'),1) = 00.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1,4')),1) = 00.0;
        end
        if(isKey(ix_,'X4,1'))
            pb.x0(ix_('X4,1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4,1')),1) = 0.0;
        end
        if(isKey(ix_,'B1'))
            pb.x0(ix_('B1'),1) = 95.0;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 95.0;
        end
        if(isKey(ix_,'X2,4'))
            pb.x0(ix_('X2,4'),1) = 0.00;
        else
            pb.y0(find(pbm.congrps==ig_('X2,4')),1) = 0.00;
        end
        if(isKey(ix_,'X4,2'))
            pb.x0(ix_('X4,2'),1) = 0.00;
        else
            pb.y0(find(pbm.congrps==ig_('X4,2')),1) = 0.00;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = 95.0;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = 95.0;
        end
        if(isKey(ix_,'X3,4'))
            pb.x0(ix_('X3,4'),1) = 0.00;
        else
            pb.y0(find(pbm.congrps==ig_('X3,4')),1) = 0.00;
        end
        if(isKey(ix_,'X4,3'))
            pb.x0(ix_('X4,3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4,3')),1) = 0.0;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = 19.0;
        end
        if(isKey(ix_,'B4'))
            pb.x0(ix_('B4'),1) = 70.0;
        else
            pb.y0(find(pbm.congrps==ig_('B4')),1) = 70.0;
        end
        if(isKey(ix_,'X5,4'))
            pb.x0(ix_('X5,4'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5,4')),1) = 0.0;
        end
        if(isKey(ix_,'X4,5'))
            pb.x0(ix_('X4,5'),1) = 00.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4,5')),1) = 00.0;
        end
        if(isKey(ix_,'X6,5'))
            pb.x0(ix_('X6,5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X6,5')),1) = 0.0;
        end
        if(isKey(ix_,'X5,6'))
            pb.x0(ix_('X5,6'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5,6')),1) = 0.0;
        end
        if(isKey(ix_,'X7,5'))
            pb.x0(ix_('X7,5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X7,5')),1) = 0.0;
        end
        if(isKey(ix_,'X5,7'))
            pb.x0(ix_('X5,7'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5,7')),1) = 0.0;
        end
        if(isKey(ix_,'B5'))
            pb.x0(ix_('B5'),1) = 70.0;
        else
            pb.y0(find(pbm.congrps==ig_('B5')),1) = 70.0;
        end
        if(isKey(ix_,'B6'))
            pb.x0(ix_('B6'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B6')),1) = 19.0;
        end
        if(isKey(ix_,'B7'))
            pb.x0(ix_('B7'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B7')),1) = 19.0;
        end
        if(isKey(ix_,'X8,5'))
            pb.x0(ix_('X8,5'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8,5')),1) = 0.0;
        end
        if(isKey(ix_,'X5,8'))
            pb.x0(ix_('X5,8'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5,8')),1) = 0.0;
        end
        if(isKey(ix_,'X9,8'))
            pb.x0(ix_('X9,8'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X9,8')),1) = 0.0;
        end
        if(isKey(ix_,'X8,9'))
            pb.x0(ix_('X8,9'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8,9')),1) = 0.0;
        end
        if(isKey(ix_,'X10,8'))
            pb.x0(ix_('X10,8'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X10,8')),1) = 0.0;
        end
        if(isKey(ix_,'X8,10'))
            pb.x0(ix_('X8,10'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8,10')),1) = 0.0;
        end
        if(isKey(ix_,'X11,8'))
            pb.x0(ix_('X11,8'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X11,8')),1) = 0.0;
        end
        if(isKey(ix_,'X8,11'))
            pb.x0(ix_('X8,11'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8,11')),1) = 0.0;
        end
        if(isKey(ix_,'B8'))
            pb.x0(ix_('B8'),1) = 70.0;
        else
            pb.y0(find(pbm.congrps==ig_('B8')),1) = 70.0;
        end
        if(isKey(ix_,'B9'))
            pb.x0(ix_('B9'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B9')),1) = 19.0;
        end
        if(isKey(ix_,'B10'))
            pb.x0(ix_('B10'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B10')),1) = 19.0;
        end
        if(isKey(ix_,'B11'))
            pb.x0(ix_('B11'),1) = 19.0;
        else
            pb.y0(find(pbm.congrps==ig_('B11')),1) = 19.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eBETA1',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eBETA2',iet_);
        elftv{it}{1} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCOMA1',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eCOMA2',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eCOMB1',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eCOMB2',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'EB1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA1';
        ielftype(ie) = iet_('eBETA1');
        vname = 'B1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA1';
        ielftype(ie) = iet_('eBETA1');
        vname = 'B2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA1';
        ielftype(ie) = iet_('eBETA1');
        vname = 'B4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA1';
        ielftype(ie) = iet_('eBETA1');
        vname = 'B5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA1';
        ielftype(ie) = iet_('eBETA1');
        vname = 'B8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EB11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eBETA2';
        ielftype(ie) = iet_('eBETA2');
        vname = 'B11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('1'):v_('P1')
            ename = ['EGA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            vname = ['X',int2str(round(v_('4C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('4C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            vname = ['X',int2str(round(v_('4C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I),',',int2str(round(v_('4C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I+3') = 3+I;
            ename = ['EGA',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            ename = ['EGA',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('4C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGA',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('4C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            ename = ['EGB',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('4C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+3')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('4C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I+6') = 6+I;
            v_('I+8') = 8+I;
            ename = ['EGA',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            ename = ['EGA',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('8C'))),',',int2str(round(v_('I+8')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGA',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+8'))),',',int2str(round(v_('8C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            ename = ['EGB',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('8C'))),',',int2str(round(v_('I+8')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+8'))),',',int2str(round(v_('8C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I+9') = 9+I;
            ename = ['EGA',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            ename = ['EGA',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+8'))),',',int2str(round(v_('8C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGA',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('8C'))),',',int2str(round(v_('I+8')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            ename = ['EGB',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+8'))),',',int2str(round(v_('8C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I+9')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('8C'))),',',int2str(round(v_('I+8')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('6C'):v_('7C')
            v_('I2') = 2*I;
            v_('I2+1') = 1+v_('I2');
            ename = ['EGA',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            ename = ['EGA',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('5C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGA',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('5C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            ename = ['EGB',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('5C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I2+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('5C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I2+2') = 2+v_('I2');
            ename = ['EGA',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMA1';
            ielftype(ie) = iet_('eCOMA1');
            ename = ['EGA',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('5C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGA',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('5C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eCOMB1';
            ielftype(ie) = iet_('eCOMB1');
            ename = ['EGB',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(I),',',int2str(round(v_('5C')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EGB',int2str(round(v_('I2+2')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('5C'))),',',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'EGA17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMA2';
        ielftype(ie) = iet_('eCOMA2');
        vname = 'X5,4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGB17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMB2';
        ielftype(ie) = iet_('eCOMB2');
        vname = 'X5,4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGA18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMA2';
        ielftype(ie) = iet_('eCOMA2');
        vname = 'X4,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5,4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGB18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMB2';
        ielftype(ie) = iet_('eCOMB2');
        vname = 'X4,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5,4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGA19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMA2';
        ielftype(ie) = iet_('eCOMA2');
        vname = 'X5,8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGB19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMB2';
        ielftype(ie) = iet_('eCOMB2');
        vname = 'X5,8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGA20';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMA2';
        ielftype(ie) = iet_('eCOMA2');
        vname = 'X8,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5,8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EGB20';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eCOMB2';
        ielftype(ie) = iet_('eCOMB2');
        vname = 'X8,5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('V',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5,8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('W',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for I=v_('1'):v_('NLINK')
            ig = ig_(['GB',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EGB',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['GA',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EGA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-MN-31-31';
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
%  SET UP THE ELEMENTS *
%  ROUTINE             *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 80.0;
        pbm.efpar(2) = 20.0;
        varargout{1} = pbm;

    case 'eBETA2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CB = 20.0;
        varargout{1} = EV_(1)/(CB-EV_(1));
        if(nargout>1)
            g_(1,1) = CB/((CB-EV_(1))^2);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2*CB/((CB-EV_(1))^3);
                varargout{3} = H_;
            end
        end

    case 'eBETA1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CB = 100.0;
        varargout{1} = EV_(1)/(CB-EV_(1));
        if(nargout>1)
            g_(1,1) = CB/((CB-EV_(1))^2);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2*CB/((CB-EV_(1))^3);
                varargout{3} = H_;
            end
        end

    case 'eCOMA1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CIJ = 1000.0;
        varargout{1} = EV_(1)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)));
        if(nargout>1)
            g_(1,1) =...
                  (CIJ-pbm.efpar(2)*EV_(2))/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^2;
            g_(2,1) =...
                  EV_(1)*pbm.efpar(2)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) =...
                      2*pbm.efpar(1)*(CIJ-pbm.efpar(2)*EV_(2))/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                H_(1,2) = pbm.efpar(2)*(CIJ+pbm.efpar(1)*EV_(1)-pbm.efpar(2)*EV_(2))/...
                     (CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                H_(2,1) = H_(1,2);
                H_(2,2) =...
                      2*pbm.efpar(2)*pbm.efpar(2)*EV_(1)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                varargout{3} = H_;
            end
        end

    case 'eCOMB1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CIJ = 1000.0;
        varargout{1} = EV_(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)));
        if(nargout>1)
            g_(1,1) =...
                  (CIJ-pbm.efpar(1)*EV_(2))/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^2;
            g_(2,1) =...
                  EV_(1)*pbm.efpar(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) =...
                      2*pbm.efpar(2)*(CIJ-pbm.efpar(1)*EV_(2))/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                H_(1,2) = pbm.efpar(1)*(CIJ+pbm.efpar(2)*EV_(1)-pbm.efpar(1)*EV_(2))/...
                     (CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                H_(2,1) = H_(1,2);
                H_(2,2) =...
                      2*pbm.efpar(1)*pbm.efpar(1)*EV_(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                varargout{3} = H_;
            end
        end

    case 'eCOMA2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CIJ = 10000.0;
        varargout{1} = EV_(1)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)));
        if(nargout>1)
            g_(1,1) =...
                  (CIJ-pbm.efpar(2)*EV_(2))/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^2;
            g_(2,1) =...
                  EV_(1)*pbm.efpar(2)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) =...
                      2*pbm.efpar(1)*(CIJ-pbm.efpar(2)*EV_(2))/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                H_(1,2) = pbm.efpar(2)*(CIJ+pbm.efpar(1)*EV_(1)-pbm.efpar(2)*EV_(2))/...
                     (CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                H_(2,1) = H_(1,2);
                H_(2,2) =...
                      2*pbm.efpar(2)*pbm.efpar(2)*EV_(1)/(CIJ-(pbm.efpar(1)*EV_(1)+pbm.efpar(2)*EV_(2)))^3;
                varargout{3} = H_;
            end
        end

    case 'eCOMB2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        CIJ = 10000.0;
        varargout{1} = EV_(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)));
        if(nargout>1)
            g_(1,1) =...
                  (CIJ-pbm.efpar(1)*EV_(2))/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^2;
            g_(2,1) =...
                  EV_(1)*pbm.efpar(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) =...
                      2*pbm.efpar(2)*(CIJ-pbm.efpar(1)*EV_(2))/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                H_(1,2) = pbm.efpar(1)*(CIJ+pbm.efpar(2)*EV_(1)-pbm.efpar(1)*EV_(2))/...
                     (CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                H_(2,1) = H_(1,2);
                H_(2,2) =...
                      2*pbm.efpar(1)*pbm.efpar(1)*EV_(1)/(CIJ-(pbm.efpar(2)*EV_(1)+pbm.efpar(1)*EV_(2)))^3;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [2,0];
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

