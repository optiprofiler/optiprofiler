function varargout = QC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : QC
%    *********
% 
%    Source: Quality Control problem 104 from
%    Betty Schultz and Ben Reiser.
% 
%    SIF input: Andrew Conn, August 1992.
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'C-COLR2-MY-9-4'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'QC';

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
        v_('M') = 9;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_('15') = 15;
        v_('16') = 16;
        v_('17') = 17;
        v_('N') = 10048;
        v_('NGG') = 9900;
        v_('NGB') = 35;
        v_('NBG') = 15;
        v_('NBU') = 18;
        v_('NBB') = 2;
        v_('NUU') = 40;
        v_('NUB') = 21;
        v_('NGBUU') = 9;
        v_('NBGUU') = 2;
        v_('NBGBU') = 2;
        v_('NGBUB') = 3;
        v_('NBBUU') = 1;
        v_('NBBBU') = 0;
        v_('NBBUB') = 0;
        v_('MM') = v_('N')-v_('NGG');
        v_('S1') = v_('NGG')+v_('NGB');
        v_('S1') = v_('S1')+v_('NGBUU');
        v_('S1') = v_('S1')+v_('NGBUB');
        v_('S2') = v_('NBG')+v_('NBB');
        v_('S2') = v_('S2')+v_('NBU');
        v_('S2') = v_('S2')+v_('NBGUU');
        v_('S2') = v_('S2')+v_('NBBUU');
        v_('S2') = v_('S2')+v_('NBGBU');
        v_('S2') = v_('S2')+v_('NBBBU');
        v_('S2') = v_('S2')+v_('NBBUB');
        v_('S3') = v_('NGG')+v_('NBG');
        v_('S3') = v_('S3')+v_('NBGUU');
        v_('S3') = v_('S3')+v_('NBGBU');
        v_('S4') = v_('NGB')+v_('NBB');
        v_('S4') = v_('S4')+v_('NUB');
        v_('S4') = v_('S4')+v_('NGBUU');
        v_('S4') = v_('S4')+v_('NBBUU');
        v_('S4') = v_('S4')+v_('NBBBU');
        v_('S4') = v_('S4')+v_('NGBUB');
        v_('S4') = v_('S4')+v_('NBBUB');
        v_('U1') = v_('NGB')+v_('NGBUU');
        v_('U1') = v_('U1')+v_('NGBUB');
        v_('L1') = v_('U1')+v_('NUU');
        v_('L1') = v_('L1')+v_('NUB');
        v_('U2') = v_('NBG')+v_('NBGUU');
        v_('U2') = v_('U2')+v_('NBGBU');
        v_('L2') = v_('U2')+v_('NUU');
        v_('L2') = v_('L2')+v_('NBU');
        v_('U3') = v_('NBB')+v_('NBBUU');
        v_('U3') = v_('U3')+v_('NBBBU');
        v_('U3') = v_('U3')+v_('NBBUB');
        v_('L3') = v_('U3')+v_('NUU');
        v_('L3') = v_('L3')+v_('NBU');
        v_('L3') = v_('L3')+v_('NUB');
        v_('U4') = v_('U2')+v_('NBU');
        v_('L4') = v_('U4')+v_('NUU');
        v_('U5') = v_('U1')+v_('NUB');
        v_('L5') = v_('U5')+v_('NUU');
        v_('U6') = v_('U3')+v_('NUB');
        v_('L6') = v_('U6')+v_('NBU');
        v_('L6') = v_('L6')+v_('NUU');
        v_('U7') = v_('U3')+v_('NBU');
        v_('L7') = v_('U7')+v_('NUB');
        v_('L7') = v_('L7')+v_('NUU');
        v_('L8') = v_('U2')+v_('NBBUU');
        v_('L8') = v_('L8')+v_('NBBBU');
        v_('L8') = v_('L8')+v_('NBBUB');
        v_('L8') = v_('L8')+v_('NBU');
        v_('L8') = v_('L8')+v_('NBB');
        v_('U8') = v_('L8')+v_('NUU');
        v_('U8') = v_('U8')+v_('NUB');
        v_('L9') = v_('U1')+v_('NBBUU');
        v_('L9') = v_('L9')+v_('NBBBU');
        v_('L9') = v_('L9')+v_('NBBUB');
        v_('L9') = v_('L9')+v_('NUB');
        v_('L9') = v_('L9')+v_('NBB');
        v_('U9') = v_('L9')+v_('NUU');
        v_('U9') = v_('U9')+v_('NBU');
        v_('ZERO') = 0.0;
        v_('TWO') = 2.0;
        v_('RS1') = v_('S1');
        v_('RS2') = v_('S2');
        v_('RS3') = v_('S3');
        v_('RS4') = v_('S4');
        v_('RNBG') = v_('NBG');
        v_('RNBGBU') = v_('NBGBU');
        v_('RNBGUU') = v_('NBGUU');
        v_('RNGB') = v_('NGB');
        v_('RNGBUB') = v_('NGBUB');
        v_('RNGBUU') = v_('NGBUU');
        v_('RNBB') = v_('NBB');
        v_('RNBBUB') = v_('NBBUB');
        v_('RNBBBU') = v_('NBBBU');
        v_('RNBBUU') = v_('NBBUU');
        v_('RNBGBU') = v_('NBGBU');
        v_('RNUU') = v_('NUU');
        v_('RNUB') = v_('NUB');
        v_('RNBU') = v_('NBU');
        v_('RN') = v_('N');
        v_('RL1') = v_('L1');
        v_('RU1') = v_('U1');
        v_('RL2') = v_('L2');
        v_('RU2') = v_('U2');
        v_('RL3') = v_('L3');
        v_('RU3') = v_('U3');
        v_('RL4') = v_('L4');
        v_('RU4') = v_('U4');
        v_('RL5') = v_('L5');
        v_('RU5') = v_('U5');
        v_('RL6') = v_('L6');
        v_('RU6') = v_('U6');
        v_('RL7') = v_('L7');
        v_('RU7') = v_('U7');
        v_('RL8') = v_('L8');
        v_('RU8') = v_('U8');
        v_('RL9') = v_('L9');
        v_('RU9') = v_('U9');
        v_('LF1') = v_('RL8')/v_('RN');
        v_('UF1') = v_('RU8')/v_('RN');
        v_('SF1') = v_('LF1')+v_('UF1');
        v_('SF1') = v_('SF1')/v_('TWO');
        v_('LF2') = v_('RL9')/v_('RN');
        v_('UF2') = v_('RU9')/v_('RN');
        v_('SF2') = v_('LF2')+v_('UF2');
        v_('SF2') = v_('SF2')/v_('TWO');
        v_('LGBGB') = v_('RNGB')/v_('RL1');
        v_('UGBGB') = v_('RNGB')/v_('RU1');
        v_('SGBGB') = v_('LGBGB')+v_('UGBGB');
        v_('SGBGB') = v_('SGBGB')/v_('TWO');
        v_('LUBGB') = v_('RNGBUB')/v_('RL5');
        v_('UUBGB') = v_('RNGBUB')+v_('RNUB');
        v_('UUBGB') = v_('UUBGB')/v_('RU5');
        v_('SUBGB') = v_('LUBGB')+v_('UUBGB');
        v_('SUBGB') = v_('SUBGB')/v_('TWO');
        v_('LBGBG') = v_('RNBG')/v_('RL2');
        v_('UBGBG') = v_('RNBG')/v_('RU2');
        v_('SBGBG') = v_('LBGBG')+v_('UBGBG');
        v_('SBGBG') = v_('SBGBG')/v_('TWO');
        v_('LBUBG') = v_('RNBGBU')/v_('RL4');
        v_('UBUBG') = v_('RNBGBU')+v_('RNBU');
        v_('UBUBG') = v_('UBUBG')/v_('RU4');
        v_('SBUBG') = v_('LBUBG')+v_('UBUBG');
        v_('SBUBG') = v_('SBUBG')/v_('TWO');
        v_('LBBBB') = v_('RNBB')/v_('RL3');
        v_('UBBBB') = v_('RNBB')/v_('RU3');
        v_('SBBBB') = v_('LBBBB')+v_('UBBBB');
        v_('SBBBB') = v_('SBBBB')/v_('TWO');
        v_('LBUBB') = v_('RNBBBU')/v_('RL7');
        v_('UBUBB') = v_('RNBBBU')+v_('RNBU');
        v_('UBUBB') = v_('UBUBB')/v_('RU7');
        v_('SBUBB') = v_('LBUBB')+v_('UBUBB');
        v_('SBUBB') = v_('SBUBB')/v_('TWO');
        v_('LUBBB') = v_('RNBBUB')/v_('RL6');
        v_('UUBBB') = v_('RNBBUB')+v_('RNUB');
        v_('UUBBB') = v_('UUBBB')/v_('RU6');
        v_('SUBBB') = v_('LUBBB')+v_('UUBBB');
        v_('SUBBB') = v_('SUBBB')/v_('TWO');
        v_('RMM') = v_('MM');
        v_('RMM') = v_('RMM')/v_('RN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','F1',ix_);
        pb.xnames{iv} = 'F1';
        [iv,ix_] = s2mpjlib('ii','F2',ix_);
        pb.xnames{iv} = 'F2';
        [iv,ix_] = s2mpjlib('ii','PBGBG',ix_);
        pb.xnames{iv} = 'PBGBG';
        [iv,ix_] = s2mpjlib('ii','PBUBG',ix_);
        pb.xnames{iv} = 'PBUBG';
        [iv,ix_] = s2mpjlib('ii','PGBGB',ix_);
        pb.xnames{iv} = 'PGBGB';
        [iv,ix_] = s2mpjlib('ii','PUBGB',ix_);
        pb.xnames{iv} = 'PUBGB';
        [iv,ix_] = s2mpjlib('ii','PBBBB',ix_);
        pb.xnames{iv} = 'PBBBB';
        [iv,ix_] = s2mpjlib('ii','PUBBB',ix_);
        pb.xnames{iv} = 'PUBBB';
        [iv,ix_] = s2mpjlib('ii','PBUBB',ix_);
        pb.xnames{iv} = 'PBUBB';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('1')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('F1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('2')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('F1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('3')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('F2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('4')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('F2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('5')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBUBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('6')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PUBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('7')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('F1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('F2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('8')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PGBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('9')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBGBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('10')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('11')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PGBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('PUBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('12')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBGBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('PBUBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('13')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('PUBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('13')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBUBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('14')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBUBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('15')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PBUBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('16')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PUBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(round(v_('17')))],ig_);
        gtype{ig} = '<>';
        iv = ix_('PUBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CON0',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CON0';
        iv = ix_('PBGBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('PBUBG');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CON1';
        iv = ix_('PGBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('PUBGB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CON2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CON2';
        iv = ix_('PBBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('PUBBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('PBUBB');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','CON3',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'CON3';
        iv = ix_('F1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('F2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        pbm.gconst(ig_(['OBJ',int2str(round(v_('1')))])) = -1.0;
        pbm.gconst(ig_(['OBJ',int2str(round(v_('3')))])) = -1.0;
        pbm.gconst(ig_(['OBJ',int2str(round(v_('7')))])) = -1.0;
        pbm.gconst(ig_(['OBJ',int2str(round(v_('11')))])) = -1.0;
        pbm.gconst(ig_(['OBJ',int2str(round(v_('12')))])) = -1.0;
        pbm.gconst(ig_(['OBJ',int2str(round(v_('13')))])) = -1.0;
        pbm.gconst(ig_('CON0')) = 1.0;
        pbm.gconst(ig_('CON1')) = 1.0;
        pbm.gconst(ig_('CON2')) = 1.0;
        pbm.gconst(ig_('CON3')) = v_('RMM');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 1.0*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('F1'),1) = v_('LF1');
        pb.xupper(ix_('F1')) = v_('UF1');
        pb.xlower(ix_('F2'),1) = v_('LF2');
        pb.xupper(ix_('F2')) = v_('UF2');
        pb.xlower(ix_('PBGBG'),1) = v_('LBGBG');
        pb.xupper(ix_('PBGBG')) = v_('UBGBG');
        pb.xlower(ix_('PBUBG'),1) = v_('LBUBG');
        pb.xupper(ix_('PBUBG')) = v_('UBUBG');
        pb.xlower(ix_('PGBGB'),1) = v_('LGBGB');
        pb.xupper(ix_('PGBGB')) = v_('UGBGB');
        pb.xlower(ix_('PUBGB'),1) = v_('LUBGB');
        pb.xupper(ix_('PUBGB')) = v_('UUBGB');
        pb.xlower(ix_('PBBBB'),1) = v_('LBBBB');
        pb.xupper(ix_('PBBBB')) = v_('UBBBB');
        pb.xlower(ix_('PBUBB'),1) = 0.0;
        pb.xupper(ix_('PBUBB'),1) = 0.0;
        pb.xlower(ix_('PUBBB'),1) = 0.0;
        pb.xupper(ix_('PUBBB'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.5*ones(pb.n,1);
        pb.x0(ix_('F1'),1) = v_('SF1');
        pb.x0(ix_('F2'),1) = v_('SF2');
        pb.x0(ix_('PBGBG'),1) = v_('SBGBG');
        pb.x0(ix_('PBUBG'),1) = v_('SBUBG');
        pb.x0(ix_('PGBGB'),1) = v_('SGBGB');
        pb.x0(ix_('PUBGB'),1) = v_('SUBGB');
        pb.x0(ix_('PBBBB'),1) = v_('SBBBB');
        pb.x0(ix_('PBUBB'),1) = v_('ZERO');
        pb.x0(ix_('PUBBB'),1) = v_('ZERO');
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        [it,iet_] = s2mpjlib( 'ii', 'eI2PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        [it,iet_] = s2mpjlib( 'ii', 'eI3PROD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        [it,iet_] = s2mpjlib( 'ii', 'en3PRODI',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PROD';
        ielftype(ie) = iet_('en2PROD');
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('1')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBUBG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PROD';
        ielftype(ie) = iet_('en2PROD');
        ename = ['E',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('2')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBUBB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PROD';
        ielftype(ie) = iet_('en2PROD');
        ename = ['E',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('3')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PUBGB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en2PROD';
        ielftype(ie) = iet_('en2PROD');
        ename = ['E',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('4')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PUBBB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eI2PROD';
        ielftype(ie) = iet_('eI2PROD');
        ename = ['E',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBGBG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('5')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBUBG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eI3PROD';
        ielftype(ie) = iet_('eI3PROD');
        ename = ['E',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBGBG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('6')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBUBG';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('7')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eI2PROD';
        ielftype(ie) = iet_('eI2PROD');
        ename = ['E',int2str(round(v_('7')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('7')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PGBGB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('7')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PUBGB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('8')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eI3PROD';
        ielftype(ie) = iet_('eI3PROD');
        ename = ['E',int2str(round(v_('8')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('8')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('8')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PGBGB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('8')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PUBGB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PRODI';
        ielftype(ie) = iet_('en3PRODI');
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'F2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBBBB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PUBBB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = ['E',int2str(round(v_('9')))];
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        vname = 'PBUBB';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],1.0,1.5);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        [it,igt_] = s2mpjlib('ii','gLOG',igt_);
        grftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_(['OBJ',int2str(round(v_('1')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RS1');
        ig = ig_(['OBJ',int2str(round(v_('2')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RS2');
        ig = ig_(['OBJ',int2str(round(v_('3')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RS3');
        ig = ig_(['OBJ',int2str(round(v_('4')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RS4');
        ig = ig_(['OBJ',int2str(round(v_('5')))]);
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBU');
        ig = ig_(['OBJ',int2str(round(v_('6')))]);
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNUB');
        ig = ig_(['OBJ',int2str(round(v_('7')))]);
        pbm.grftype{ig} = 'gLOG';
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('6')))]);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('7')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('8')))]);
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('9')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNUU');
        ig = ig_(['OBJ',int2str(round(v_('8')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNGB');
        ig = ig_(['OBJ',int2str(round(v_('9')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBG');
        ig = ig_(['OBJ',int2str(round(v_('10')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBB');
        ig = ig_(['OBJ',int2str(round(v_('11')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNGBUU');
        ig = ig_(['OBJ',int2str(round(v_('12')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBGUU');
        ig = ig_(['OBJ',int2str(round(v_('13')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBBUU');
        ig = ig_(['OBJ',int2str(round(v_('14')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBGBU');
        ig = ig_(['OBJ',int2str(round(v_('15')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBBBU');
        ig = ig_(['OBJ',int2str(round(v_('16')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNGBUB');
        ig = ig_(['OBJ',int2str(round(v_('17')))]);
        pbm.grftype{ig} = 'gLOG';
        [~,posgp] = ismember('P',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('RNBBUB');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               1138.416240
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-MY-9-4';
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

    case 'en2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 0.0;
                H_(2,2) = 0.0;
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eI2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)-1;
        U_(2,3) = U_(2,3)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 0.0;
                H_(2,2) = 0.0;
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eI3PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(3,3) = U_(3,3)-1;
        U_(3,4) = U_(3,4)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*(1.0+IV_(3));
        if(nargout>1)
            g_(1,1) = IV_(2)*(1.0+IV_(3));
            g_(2,1) = IV_(1)*(1.0+IV_(3));
            g_(3,1) = IV_(1)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 0.0;
                H_(1,2) = 1.0+IV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = 0.0;
                H_(2,3) = IV_(1);
                H_(3,2) = H_(2,3);
                H_(3,3) = 0.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'en3PRODI'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,5);
        U_(1,1) = U_(1,1)+1;
        U_(2,2) = U_(2,2)+1;
        U_(3,3) = U_(3,3)-1;
        U_(3,4) = U_(3,4)-1;
        U_(3,5) = U_(3,5)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        varargout{1} = IV_(1)*IV_(2)*(1.0+IV_(3));
        if(nargout>1)
            g_(1,1) = IV_(2)*(1.0+IV_(3));
            g_(2,1) = IV_(1)*(1.0+IV_(3));
            g_(3,1) = IV_(1)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 0.0;
                H_(1,2) = 1.0+IV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = 0.0;
                H_(2,3) = IV_(1);
                H_(3,2) = H_(2,3);
                H_(3,3) = 0.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gLOG'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        T = abs(GVAR_);
        SMALL = 1.0e-10;
        LARGE = 1.0e+10;
        ARG0 = T<=SMALL;
        if(ARG0)
            FF = pbm.grpar{igr_}(1)*log(SMALL);
        end
        if(~ARG0)
            FF = pbm.grpar{igr_}(1)*log(T);
        end
        if(ARG0)
            GG = pbm.grpar{igr_}(1)*LARGE;
        end
        if(~ARG0)
            GG = pbm.grpar{igr_}(1)/T;
        end
        if(ARG0)
            HH = -pbm.grpar{igr_}(1)*LARGE^2;
        end
        if(~ARG0)
            HH = -pbm.grpar{igr_}(1)/T^2;
        end
        varargout{1} = FF;
        if(nargout>1)
            g_ = GG;
            varargout{2} = g_;
            if(nargout>2)
                H_ = HH;
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

