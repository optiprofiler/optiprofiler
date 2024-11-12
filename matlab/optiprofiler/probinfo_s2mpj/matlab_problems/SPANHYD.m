function varargout = SPANHYD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPANHYD
%    *********
% 
%    The spanish hydro-electric reservoir management problem
%    The problem is also named "HYD33" in some references.
%    This is a nonlinear network problem
% 
%    Source:
%    J.L. de la Fuente,
%    "Programacion no-lineal: applicationes en analisis, gestion y
%    planificacion de sistemas electricos",
%    Hidroelectrica Espanola, private communication, 1986.
% 
%    SIF input: Ph. Toint, Sept 1990.
% 
%    classification = 'C-CONR2-RN-97-33'
% 
%    Number of arcs    = 97
%    Number of nodes
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SPANHYD';

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
        v_('NODES') = 33;
        v_('1') = 1;
        v_('P1A') = 85.459510678;
        v_('P2A') = -6.255270637;
        v_('P3A') = 1.0862222572;
        v_('P1B') = 81.978824852;
        v_('P2B') = -6.021899007;
        v_('P3B') = 8.4278676858;
        v_('P1C') = 83.420024838;
        v_('P2C') = -6.089102950;
        v_('P3C') = 8.9209611766;
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('NODES')
            [ig,ig_] = s2mpjlib('ii',['N',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['N',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = ig_.Count;
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X5',ix_);
        pb.xnames{iv} = 'X5';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X6',ix_);
        pb.xnames{iv} = 'X6';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X7',ix_);
        pb.xnames{iv} = 'X7';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X8',ix_);
        pb.xnames{iv} = 'X8';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X9',ix_);
        pb.xnames{iv} = 'X9';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X10',ix_);
        pb.xnames{iv} = 'X10';
        ig = ig_('N32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X11',ix_);
        pb.xnames{iv} = 'X11';
        ig = ig_('N1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X12',ix_);
        pb.xnames{iv} = 'X12';
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X13',ix_);
        pb.xnames{iv} = 'X13';
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X14',ix_);
        pb.xnames{iv} = 'X14';
        ig = ig_('N4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X15',ix_);
        pb.xnames{iv} = 'X15';
        ig = ig_('N5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X16',ix_);
        pb.xnames{iv} = 'X16';
        ig = ig_('N6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X17',ix_);
        pb.xnames{iv} = 'X17';
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X18',ix_);
        pb.xnames{iv} = 'X18';
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X19',ix_);
        pb.xnames{iv} = 'X19';
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X20',ix_);
        pb.xnames{iv} = 'X20';
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X21',ix_);
        pb.xnames{iv} = 'X21';
        ig = ig_('N1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X22',ix_);
        pb.xnames{iv} = 'X22';
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X23',ix_);
        pb.xnames{iv} = 'X23';
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X24',ix_);
        pb.xnames{iv} = 'X24';
        ig = ig_('N4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X25',ix_);
        pb.xnames{iv} = 'X25';
        ig = ig_('N6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X26',ix_);
        pb.xnames{iv} = 'X26';
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X27',ix_);
        pb.xnames{iv} = 'X27';
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X28',ix_);
        pb.xnames{iv} = 'X28';
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X29',ix_);
        pb.xnames{iv} = 'X29';
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X30',ix_);
        pb.xnames{iv} = 'X30';
        ig = ig_('N1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X31',ix_);
        pb.xnames{iv} = 'X31';
        ig = ig_('N2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X32',ix_);
        pb.xnames{iv} = 'X32';
        ig = ig_('N3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X33',ix_);
        pb.xnames{iv} = 'X33';
        ig = ig_('N4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X34',ix_);
        pb.xnames{iv} = 'X34';
        ig = ig_('N5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X35',ix_);
        pb.xnames{iv} = 'X35';
        ig = ig_('N6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X36',ix_);
        pb.xnames{iv} = 'X36';
        ig = ig_('N7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X37',ix_);
        pb.xnames{iv} = 'X37';
        ig = ig_('N8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X38',ix_);
        pb.xnames{iv} = 'X38';
        ig = ig_('N9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X39',ix_);
        pb.xnames{iv} = 'X39';
        ig = ig_('N10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X40',ix_);
        pb.xnames{iv} = 'X40';
        ig = ig_('N11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X41',ix_);
        pb.xnames{iv} = 'X41';
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X42',ix_);
        pb.xnames{iv} = 'X42';
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X43',ix_);
        pb.xnames{iv} = 'X43';
        ig = ig_('N14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X44',ix_);
        pb.xnames{iv} = 'X44';
        ig = ig_('N15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X45',ix_);
        pb.xnames{iv} = 'X45';
        ig = ig_('N16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X46',ix_);
        pb.xnames{iv} = 'X46';
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X47',ix_);
        pb.xnames{iv} = 'X47';
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X48',ix_);
        pb.xnames{iv} = 'X48';
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X49',ix_);
        pb.xnames{iv} = 'X49';
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X50',ix_);
        pb.xnames{iv} = 'X50';
        ig = ig_('N11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X51',ix_);
        pb.xnames{iv} = 'X51';
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X52',ix_);
        pb.xnames{iv} = 'X52';
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X53',ix_);
        pb.xnames{iv} = 'X53';
        ig = ig_('N14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X54',ix_);
        pb.xnames{iv} = 'X54';
        ig = ig_('N16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X55',ix_);
        pb.xnames{iv} = 'X55';
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X56',ix_);
        pb.xnames{iv} = 'X56';
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X57',ix_);
        pb.xnames{iv} = 'X57';
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X58',ix_);
        pb.xnames{iv} = 'X58';
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X59',ix_);
        pb.xnames{iv} = 'X59';
        ig = ig_('N11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X60',ix_);
        pb.xnames{iv} = 'X60';
        ig = ig_('N12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X61',ix_);
        pb.xnames{iv} = 'X61';
        ig = ig_('N13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X62',ix_);
        pb.xnames{iv} = 'X62';
        ig = ig_('N14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X63',ix_);
        pb.xnames{iv} = 'X63';
        ig = ig_('N15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X64',ix_);
        pb.xnames{iv} = 'X64';
        ig = ig_('N16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X65',ix_);
        pb.xnames{iv} = 'X65';
        ig = ig_('N17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X66',ix_);
        pb.xnames{iv} = 'X66';
        ig = ig_('N18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X67',ix_);
        pb.xnames{iv} = 'X67';
        ig = ig_('N19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X68',ix_);
        pb.xnames{iv} = 'X68';
        ig = ig_('N20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X69',ix_);
        pb.xnames{iv} = 'X69';
        ig = ig_('N21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X70',ix_);
        pb.xnames{iv} = 'X70';
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X71',ix_);
        pb.xnames{iv} = 'X71';
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X72',ix_);
        pb.xnames{iv} = 'X72';
        ig = ig_('N24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X73',ix_);
        pb.xnames{iv} = 'X73';
        ig = ig_('N25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X74',ix_);
        pb.xnames{iv} = 'X74';
        ig = ig_('N26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X75',ix_);
        pb.xnames{iv} = 'X75';
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X76',ix_);
        pb.xnames{iv} = 'X76';
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X77',ix_);
        pb.xnames{iv} = 'X77';
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X78',ix_);
        pb.xnames{iv} = 'X78';
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X79',ix_);
        pb.xnames{iv} = 'X79';
        ig = ig_('N21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X80',ix_);
        pb.xnames{iv} = 'X80';
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X81',ix_);
        pb.xnames{iv} = 'X81';
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X82',ix_);
        pb.xnames{iv} = 'X82';
        ig = ig_('N24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X83',ix_);
        pb.xnames{iv} = 'X83';
        ig = ig_('N26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X84',ix_);
        pb.xnames{iv} = 'X84';
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X85',ix_);
        pb.xnames{iv} = 'X85';
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X86',ix_);
        pb.xnames{iv} = 'X86';
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X87',ix_);
        pb.xnames{iv} = 'X87';
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X88',ix_);
        pb.xnames{iv} = 'X88';
        ig = ig_('N21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X89',ix_);
        pb.xnames{iv} = 'X89';
        ig = ig_('N22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X90',ix_);
        pb.xnames{iv} = 'X90';
        ig = ig_('N23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X91',ix_);
        pb.xnames{iv} = 'X91';
        ig = ig_('N24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X92',ix_);
        pb.xnames{iv} = 'X92';
        ig = ig_('N25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X93',ix_);
        pb.xnames{iv} = 'X93';
        ig = ig_('N26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X94',ix_);
        pb.xnames{iv} = 'X94';
        ig = ig_('N27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X95',ix_);
        pb.xnames{iv} = 'X95';
        ig = ig_('N28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X96',ix_);
        pb.xnames{iv} = 'X96';
        ig = ig_('N29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X97',ix_);
        pb.xnames{iv} = 'X97';
        ig = ig_('N30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('N33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        pbm.gconst(ig_('N1')) = -5.13800e+01;
        pbm.gconst(ig_('N2')) = -1.38400e+01;
        pbm.gconst(ig_('N3')) = -2.58000;
        pbm.gconst(ig_('N4')) = -2.19100e+01;
        pbm.gconst(ig_('N6')) = -1.29700e+01;
        pbm.gconst(ig_('N8')) = -2.89000;
        pbm.gconst(ig_('N9')) = -2.08400e+01;
        pbm.gconst(ig_('N10')) = -1.71400e+01;
        pbm.gconst(ig_('N11')) = -3.20600e+01;
        pbm.gconst(ig_('N12')) = -2.80000e-01;
        pbm.gconst(ig_('N13')) = -4.20000;
        pbm.gconst(ig_('N14')) = -4.83700e+01;
        pbm.gconst(ig_('N16')) = -1.81300e+01;
        pbm.gconst(ig_('N18')) = 1.61000;
        pbm.gconst(ig_('N19')) = -2.66000e+01;
        pbm.gconst(ig_('N20')) = -1.87600e+01;
        pbm.gconst(ig_('N21')) = -1.81300e+01;
        pbm.gconst(ig_('N24')) = -1.81300e+01;
        pbm.gconst(ig_('N26')) = -9.10000;
        pbm.gconst(ig_('N28')) = 5.81000;
        pbm.gconst(ig_('N29')) = -9.10000;
        pbm.gconst(ig_('N30')) = -6.02000;
        pbm.gconst(ig_('N31')) = 6.08350e+02;
        pbm.gconst(ig_('N32')) = -4.62634e+03;
        pbm.gconst(ig_('N33')) = 4.36300e+03;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 3.02400e+03*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('X1'),1) = 7.70000e+01;
        pb.xupper(ix_('X1')) = 7.70100e+01;
        pb.xlower(ix_('X2'),1) = 1.12452e+03;
        pb.xupper(ix_('X2')) = 1.12453e+03;
        pb.xlower(ix_('X3'),1) = 1.58000e+02;
        pb.xupper(ix_('X3')) = 1.58010e+02;
        pb.xlower(ix_('X4'),1) = 1.60000e+01;
        pb.xupper(ix_('X4')) = 1.60100e+01;
        pb.xlower(ix_('X5'),1) = 0.00000;
        pb.xupper(ix_('X5'),1) = 0.00000;
        pb.xlower(ix_('X6'),1) = 7.83650e+02;
        pb.xupper(ix_('X6')) = 7.83660e+02;
        pb.xlower(ix_('X7'),1) = 1.10000e+01;
        pb.xupper(ix_('X7')) = 1.10100e+01;
        pb.xlower(ix_('X8'),1) = 4.90000e+01;
        pb.xupper(ix_('X8')) = 4.90100e+01;
        pb.xlower(ix_('X9'),1) = 2.15517e+03;
        pb.xupper(ix_('X9')) = 2.15518e+03;
        pb.xlower(ix_('X10'),1) = 2.52000e+02;
        pb.xupper(ix_('X10')) = 2.52010e+02;
        pb.xupper(ix_('X11')) = 3.97840e+02;
        pb.xupper(ix_('X12')) = 2.22320e+02;
        pb.xupper(ix_('X13')) = 2.05630e+02;
        pb.xupper(ix_('X14')) = 2.05630e+02;
        pb.xupper(ix_('X15')) = 2.05630e+02;
        pb.xupper(ix_('X16')) = 1.24830e+02;
        pb.xupper(ix_('X17')) = 1.27010e+02;
        pb.xupper(ix_('X18')) = 6.10800e+01;
        pb.xupper(ix_('X19')) = 6.14840e+02;
        pb.xupper(ix_('X20')) = 7.78080e+02;
        pb.xupper(ix_('X25')) = 7.25760e+03;
        pb.xupper(ix_('X26')) = 1.20960e+03;
        pb.xupper(ix_('X27')) = 9.07200e+02;
        pb.xupper(ix_('X28')) = 7.25760e+03;
        pb.xupper(ix_('X29')) = 7.25760e+03;
        pb.xlower(ix_('X30'),1) = 7.70000e+01;
        pb.xupper(ix_('X30'),1) = 7.70000e+01;
        pb.xlower(ix_('X31'),1) = 4.03400e+02;
        pb.xupper(ix_('X31')) = 1.31200e+03;
        pb.xlower(ix_('X32'),1) = 1.58000e+02;
        pb.xupper(ix_('X32'),1) = 1.58000e+02;
        pb.xlower(ix_('X33'),1) = 1.60000e+01;
        pb.xupper(ix_('X33'),1) = 1.60000e+01;
        pb.xlower(ix_('X34'),1) = 0.00000;
        pb.xupper(ix_('X34'),1) = 0.00000;
        pb.xlower(ix_('X35'),1) = 5.02000e+02;
        pb.xupper(ix_('X35')) = 9.28460e+02;
        pb.xlower(ix_('X36'),1) = 1.10000e+01;
        pb.xupper(ix_('X36'),1) = 1.10000e+01;
        pb.xlower(ix_('X37'),1) = 4.90000e+01;
        pb.xupper(ix_('X37'),1) = 4.90000e+01;
        pb.xlower(ix_('X38'),1) = 9.15300e+02;
        pb.xupper(ix_('X38')) = 2.61160e+03;
        pb.xlower(ix_('X39'),1) = 2.52000e+02;
        pb.xupper(ix_('X39'),1) = 2.52000e+02;
        pb.xupper(ix_('X40')) = 3.97840e+02;
        pb.xupper(ix_('X41')) = 2.22320e+02;
        pb.xupper(ix_('X42')) = 2.05630e+02;
        pb.xupper(ix_('X43')) = 2.05630e+02;
        pb.xupper(ix_('X44')) = 2.05630e+02;
        pb.xupper(ix_('X45')) = 1.24830e+02;
        pb.xupper(ix_('X46')) = 1.27010e+02;
        pb.xupper(ix_('X47')) = 6.10800e+01;
        pb.xupper(ix_('X48')) = 6.14840e+02;
        pb.xupper(ix_('X49')) = 7.78080e+02;
        pb.xupper(ix_('X54')) = 7.25760e+03;
        pb.xupper(ix_('X55')) = 1.20960e+03;
        pb.xupper(ix_('X56')) = 9.07200e+02;
        pb.xupper(ix_('X57')) = 7.25760e+03;
        pb.xupper(ix_('X58')) = 7.25760e+03;
        pb.xlower(ix_('X59'),1) = 7.70000e+01;
        pb.xupper(ix_('X59'),1) = 7.70000e+01;
        pb.xlower(ix_('X60'),1) = 4.03400e+02;
        pb.xupper(ix_('X60')) = 1.31200e+03;
        pb.xlower(ix_('X61'),1) = 1.58000e+02;
        pb.xupper(ix_('X61'),1) = 1.58000e+02;
        pb.xlower(ix_('X62'),1) = 1.60000e+01;
        pb.xupper(ix_('X62'),1) = 1.60000e+01;
        pb.xlower(ix_('X63'),1) = 0.00000;
        pb.xupper(ix_('X63'),1) = 0.00000;
        pb.xlower(ix_('X64'),1) = 5.05640e+02;
        pb.xupper(ix_('X64')) = 9.28460e+02;
        pb.xlower(ix_('X65'),1) = 1.10000e+01;
        pb.xupper(ix_('X65'),1) = 1.10000e+01;
        pb.xlower(ix_('X66'),1) = 4.90000e+01;
        pb.xupper(ix_('X66'),1) = 4.90000e+01;
        pb.xlower(ix_('X67'),1) = 9.15300e+02;
        pb.xupper(ix_('X67')) = 2.61160e+03;
        pb.xlower(ix_('X68'),1) = 2.52000e+02;
        pb.xupper(ix_('X68'),1) = 2.52000e+02;
        pb.xupper(ix_('X69')) = 3.97840e+02;
        pb.xupper(ix_('X70')) = 2.22320e+02;
        pb.xupper(ix_('X71')) = 2.05630e+02;
        pb.xupper(ix_('X72')) = 2.05630e+02;
        pb.xupper(ix_('X73')) = 2.05630e+02;
        pb.xupper(ix_('X74')) = 1.24830e+02;
        pb.xupper(ix_('X75')) = 1.27010e+02;
        pb.xupper(ix_('X76')) = 6.10800e+01;
        pb.xupper(ix_('X77')) = 6.14840e+02;
        pb.xupper(ix_('X78')) = 7.78080e+02;
        pb.xupper(ix_('X83')) = 7.25760e+03;
        pb.xupper(ix_('X84')) = 1.20960e+03;
        pb.xupper(ix_('X85')) = 9.07200e+02;
        pb.xupper(ix_('X86')) = 7.25760e+03;
        pb.xupper(ix_('X87')) = 7.25760e+03;
        pb.xlower(ix_('X88'),1) = 7.70000e+01;
        pb.xupper(ix_('X88')) = 7.70100e+01;
        pb.xlower(ix_('X89'),1) = 1.10000e+03;
        pb.xupper(ix_('X89')) = 1.10001e+03;
        pb.xlower(ix_('X90'),1) = 1.58000e+02;
        pb.xupper(ix_('X90')) = 1.58010e+02;
        pb.xlower(ix_('X91'),1) = 1.60000e+01;
        pb.xupper(ix_('X91')) = 1.60100e+01;
        pb.xlower(ix_('X92'),1) = 0.00000;
        pb.xupper(ix_('X92'),1) = 0.00000;
        pb.xlower(ix_('X93'),1) = 7.00000e+02;
        pb.xupper(ix_('X93')) = 7.00010e+02;
        pb.xlower(ix_('X94'),1) = 1.10000e+01;
        pb.xupper(ix_('X94')) = 1.10100e+01;
        pb.xlower(ix_('X95'),1) = 4.90000e+01;
        pb.xupper(ix_('X95')) = 4.90100e+01;
        pb.xlower(ix_('X96'),1) = 2.00000e+03;
        pb.xupper(ix_('X96')) = 2.00001e+03;
        pb.xlower(ix_('X97'),1) = 2.52000e+02;
        pb.xupper(ix_('X97')) = 2.52010e+02;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 7.70000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 7.70000e+01;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 1.12452e+03;
        else
            pb.y0(find(pbm.congrps==ig('X2')),1) = 1.12452e+03;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 1.58000e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 1.58000e+02;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 1.60000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X4')),1) = 1.60000e+01;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 7.83650e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 7.83650e+02;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 1.10000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X7')),1) = 1.10000e+01;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 4.90000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X8')),1) = 4.90000e+01;
        end
        if(isKey(ix_,'X9'))
            pb.x0(ix_('X9'),1) = 2.15517e+03;
        else
            pb.y0(find(pbm.congrps==ig('X9')),1) = 2.15517e+03;
        end
        if(isKey(ix_,'X10'))
            pb.x0(ix_('X10'),1) = 2.52000e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X10')),1) = 2.52000e+02;
        end
        if(isKey(ix_,'X11'))
            pb.x0(ix_('X11'),1) = 5.13800e+01;
        else
            pb.y0(find(pbm.congrps==ig('X11')),1) = 5.13800e+01;
        end
        if(isKey(ix_,'X12'))
            pb.x0(ix_('X12'),1) = 1.40210e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X12')),1) = 1.40210e+02;
        end
        if(isKey(ix_,'X13'))
            pb.x0(ix_('X13'),1) = 1.42790e+02;
        else
            pb.y0(find(pbm.congrps==ig('X13')),1) = 1.42790e+02;
        end
        if(isKey(ix_,'X14'))
            pb.x0(ix_('X14'),1) = 2.19100e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X14')),1) = 2.19100e+01;
        end
        if(isKey(ix_,'X15'))
            pb.x0(ix_('X15'),1) = 1.64700e+02;
        else
            pb.y0(find(pbm.congrps==ig('X15')),1) = 1.64700e+02;
        end
        if(isKey(ix_,'X16'))
            pb.x0(ix_('X16'),1) = 5.81900e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X16')),1) = 5.81900e+01;
        end
        if(isKey(ix_,'X17'))
            pb.x0(ix_('X17'),1) = 5.81900e+01;
        else
            pb.y0(find(pbm.congrps==ig('X17')),1) = 5.81900e+01;
        end
        if(isKey(ix_,'X18'))
            pb.x0(ix_('X18'),1) = 6.10800e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X18')),1) = 6.10800e+01;
        end
        if(isKey(ix_,'X19'))
            pb.x0(ix_('X19'),1) = 5.66430e+02;
        else
            pb.y0(find(pbm.congrps==ig('X19')),1) = 5.66430e+02;
        end
        if(isKey(ix_,'X20'))
            pb.x0(ix_('X20'),1) = 5.83570e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X20')),1) = 5.83570e+02;
        end
        if(isKey(ix_,'X30'))
            pb.x0(ix_('X30'),1) = 7.70000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X30')),1) = 7.70000e+01;
        end
        if(isKey(ix_,'X31'))
            pb.x0(ix_('X31'),1) = 1.04953e+03;
        else
            pb.y0(find(pbm.congrps==ig_('X31')),1) = 1.04953e+03;
        end
        if(isKey(ix_,'X32'))
            pb.x0(ix_('X32'),1) = 1.58000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X32')),1) = 1.58000e+02;
        end
        if(isKey(ix_,'X33'))
            pb.x0(ix_('X33'),1) = 1.60000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X33')),1) = 1.60000e+01;
        end
        if(isKey(ix_,'X35'))
            pb.x0(ix_('X35'),1) = 7.38430e+02;
        else
            pb.y0(find(pbm.congrps==ig('X35')),1) = 7.38430e+02;
        end
        if(isKey(ix_,'X36'))
            pb.x0(ix_('X36'),1) = 1.10000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X36')),1) = 1.10000e+01;
        end
        if(isKey(ix_,'X37'))
            pb.x0(ix_('X37'),1) = 4.90000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X37')),1) = 4.90000e+01;
        end
        if(isKey(ix_,'X38'))
            pb.x0(ix_('X38'),1) = 1.83536e+03;
        else
            pb.y0(find(pbm.congrps==ig_('X38')),1) = 1.83536e+03;
        end
        if(isKey(ix_,'X39'))
            pb.x0(ix_('X39'),1) = 2.52000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X39')),1) = 2.52000e+02;
        end
        if(isKey(ix_,'X40'))
            pb.x0(ix_('X40'),1) = 3.20600e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X40')),1) = 3.20600e+01;
        end
        if(isKey(ix_,'X42'))
            pb.x0(ix_('X42'),1) = 4.20000;
        else
            pb.y0(find(pbm.congrps==ig('X42')),1) = 4.20000;
        end
        if(isKey(ix_,'X43'))
            pb.x0(ix_('X43'),1) = 4.83700e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X43')),1) = 4.83700e+01;
        end
        if(isKey(ix_,'X44'))
            pb.x0(ix_('X44'),1) = 5.25700e+01;
        else
            pb.y0(find(pbm.congrps==ig('X44')),1) = 5.25700e+01;
        end
        if(isKey(ix_,'X45'))
            pb.x0(ix_('X45'),1) = 5.98500e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X45')),1) = 5.98500e+01;
        end
        if(isKey(ix_,'X46'))
            pb.x0(ix_('X46'),1) = 5.98500e+01;
        else
            pb.y0(find(pbm.congrps==ig('X46')),1) = 5.98500e+01;
        end
        if(isKey(ix_,'X47'))
            pb.x0(ix_('X47'),1) = 5.82400e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X47')),1) = 5.82400e+01;
        end
        if(isKey(ix_,'X49'))
            pb.x0(ix_('X49'),1) = 1.87600e+01;
        else
            pb.y0(find(pbm.congrps==ig('X49')),1) = 1.87600e+01;
        end
        if(isKey(ix_,'X59'))
            pb.x0(ix_('X59'),1) = 7.70000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X59')),1) = 7.70000e+01;
        end
        if(isKey(ix_,'X60'))
            pb.x0(ix_('X60'),1) = 1.08187e+03;
        else
            pb.y0(find(pbm.congrps==ig('X60')),1) = 1.08187e+03;
        end
        if(isKey(ix_,'X61'))
            pb.x0(ix_('X61'),1) = 1.58000e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X61')),1) = 1.58000e+02;
        end
        if(isKey(ix_,'X62'))
            pb.x0(ix_('X62'),1) = 1.60000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X62')),1) = 1.60000e+01;
        end
        if(isKey(ix_,'X64'))
            pb.x0(ix_('X64'),1) = 6.96710e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X64')),1) = 6.96710e+02;
        end
        if(isKey(ix_,'X65'))
            pb.x0(ix_('X65'),1) = 1.10000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X65')),1) = 1.10000e+01;
        end
        if(isKey(ix_,'X66'))
            pb.x0(ix_('X66'),1) = 4.90000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X66')),1) = 4.90000e+01;
        end
        if(isKey(ix_,'X67'))
            pb.x0(ix_('X67'),1) = 1.97277e+03;
        else
            pb.y0(find(pbm.congrps==ig('X67')),1) = 1.97277e+03;
        end
        if(isKey(ix_,'X68'))
            pb.x0(ix_('X68'),1) = 2.52000e+02;
        else
            pb.y0(find(pbm.congrps==ig_('X68')),1) = 2.52000e+02;
        end
        if(isKey(ix_,'X69'))
            pb.x0(ix_('X69'),1) = 1.81300e+01;
        else
            pb.y0(find(pbm.congrps==ig('X69')),1) = 1.81300e+01;
        end
        if(isKey(ix_,'X72'))
            pb.x0(ix_('X72'),1) = 1.81300e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X72')),1) = 1.81300e+01;
        end
        if(isKey(ix_,'X73'))
            pb.x0(ix_('X73'),1) = 1.81300e+01;
        else
            pb.y0(find(pbm.congrps==ig('X73')),1) = 1.81300e+01;
        end
        if(isKey(ix_,'X74'))
            pb.x0(ix_('X74'),1) = 5.81000;
        else
            pb.y0(find(pbm.congrps==ig_('X74')),1) = 5.81000;
        end
        if(isKey(ix_,'X75'))
            pb.x0(ix_('X75'),1) = 5.81000;
        else
            pb.y0(find(pbm.congrps==ig('X75')),1) = 5.81000;
        end
        if(isKey(ix_,'X78'))
            pb.x0(ix_('X78'),1) = 6.02000;
        else
            pb.y0(find(pbm.congrps==ig_('X78')),1) = 6.02000;
        end
        if(isKey(ix_,'X88'))
            pb.x0(ix_('X88'),1) = 7.70000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X88')),1) = 7.70000e+01;
        end
        if(isKey(ix_,'X89'))
            pb.x0(ix_('X89'),1) = 1.10000e+03;
        else
            pb.y0(find(pbm.congrps==ig_('X89')),1) = 1.10000e+03;
        end
        if(isKey(ix_,'X90'))
            pb.x0(ix_('X90'),1) = 1.58000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X90')),1) = 1.58000e+02;
        end
        if(isKey(ix_,'X91'))
            pb.x0(ix_('X91'),1) = 1.60000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X91')),1) = 1.60000e+01;
        end
        if(isKey(ix_,'X93'))
            pb.x0(ix_('X93'),1) = 7.00000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X93')),1) = 7.00000e+02;
        end
        if(isKey(ix_,'X94'))
            pb.x0(ix_('X94'),1) = 1.10000e+01;
        else
            pb.y0(find(pbm.congrps==ig_('X94')),1) = 1.10000e+01;
        end
        if(isKey(ix_,'X95'))
            pb.x0(ix_('X95'),1) = 4.90000e+01;
        else
            pb.y0(find(pbm.congrps==ig('X95')),1) = 4.90000e+01;
        end
        if(isKey(ix_,'X96'))
            pb.x0(ix_('X96'),1) = 2.00000e+03;
        else
            pb.y0(find(pbm.congrps==ig_('X96')),1) = 2.00000e+03;
        end
        if(isKey(ix_,'X97'))
            pb.x0(ix_('X97'),1) = 2.52000e+02;
        else
            pb.y0(find(pbm.congrps==ig('X97')),1) = 2.52000e+02;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSPHYD',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V6';
        elftv{it}{6} = 'V16';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        elftv{it}{9} = 'V9';
        elftv{it}{10} = 'V10';
        elftv{it}{11} = 'V11';
        elftv{it}{12} = 'V12';
        elftv{it}{13} = 'V13';
        elftv{it}{14} = 'V14';
        elftv{it}{15} = 'V17';
        elftv{it}{16} = 'V18';
        elftv{it}{17} = 'V19';
        elftv{it}{18} = 'V20';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'EA';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSPHYD';
            ielftype(ie) = iet_('eSPHYD');
        end
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V11',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V12',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V13',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V14',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V16',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V17',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V18',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V19',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V20',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P1A');
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P2A');
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P3A');
        ename = 'EB';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSPHYD';
            ielftype(ie) = iet_('eSPHYD');
        end
        vname = 'X30';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X31';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X33';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X35';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X37';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X38';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X39';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X40';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V11',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V12',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X42';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V13',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X43';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V14',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X45';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V16',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X46';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V17',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X47';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V18',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X48';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V19',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V20',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P1B');
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P2B');
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P3B');
        ename = 'EC';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSPHYD';
            ielftype(ie) = iet_('eSPHYD');
        end
        vname = 'X59';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X60';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X61';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X62';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X64';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X65';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X66';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X67';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X68';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X69';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V11',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X70';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V12',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X71';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V13',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X72';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V14',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X74';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V16',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X75';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V17',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X76';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V18',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X77';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V19',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X78';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],3.02400e+03,0.0);
        posev = find(strcmp('V20',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P1C');
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P2C');
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('P3C');
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EA');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('EB');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EC');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               239.738001
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CONR2-RN-97-33';
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
        pbm.efpar(1) = -637.993;
        pbm.efpar(2) = -15.4452;
        pbm.efpar(3) = -6.2597;
        pbm.efpar(4) = -2.8699;
        pbm.efpar(5) = -8.6773;
        pbm.efpar(6) = -8.9261;
        pbm.efpar(7) = -672.26;
        pbm.efpar(8) = -334.117;
        pbm.efpar(9) = -205.321;
        pbm.efpar(10) = 2.02;
        pbm.efpar(11) = 2.57;
        pbm.efpar(12) = 2.2146;
        pbm.efpar(13) = 2.5128;
        pbm.efpar(14) = 2.5707;
        pbm.efpar(15) = 2.7601;
        pbm.efpar(16) = 2.5402;
        pbm.efpar(17) = 2.711;
        pbm.efpar(18) = 2.6673;
        pbm.efpar(19) = 344.904;
        pbm.efpar(20) = 283.67;
        pbm.efpar(21) = 188.597;
        pbm.efpar(22) = 212.864;
        pbm.efpar(23) = 353.204;
        pbm.efpar(24) = 316.895;
        pbm.efpar(25) = 295.055;
        pbm.efpar(26) = 163.013;
        pbm.efpar(27) = 107.338;
        pbm.efpar(28) = 0.1265;
        pbm.efpar(29) = 0.031;
        pbm.efpar(30) = 0.558;
        pbm.efpar(31) = 0.7584;
        pbm.efpar(32) = 0.0541;
        pbm.efpar(33) = 0.9038;
        pbm.efpar(34) = 0.1557;
        pbm.efpar(35) = 0.0262;
        pbm.efpar(36) = -0.3111;
        pbm.efpar(37) = -2.5225e-4;
        pbm.efpar(38) = -7.2e-6;
        pbm.efpar(39) = -1.4052e-3;
        pbm.efpar(40) = -1.4353e-2;
        pbm.efpar(41) = -2.01e-5;
        pbm.efpar(42) = 1.4139e-3;
        pbm.efpar(43) = 1.7195e-3;
        pbm.efpar(44) = -2.4e-6;
        pbm.efpar(45) = 2.372e-4;
        varargout{1} = pbm;

    case 'eSPHYD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        Z1 = pbm.efpar(19)+pbm.efpar(28)*EV_(1)+pbm.efpar(37)*EV_(1)*EV_(1);
        Z2 = pbm.efpar(20)+pbm.efpar(29)*EV_(2)+pbm.efpar(38)*EV_(2)*EV_(2);
        Z3 = pbm.efpar(21)+pbm.efpar(30)*EV_(3)+pbm.efpar(39)*EV_(3)*EV_(3);
        Z4 = pbm.efpar(22)+pbm.efpar(31)*EV_(4)+pbm.efpar(40)*EV_(4)*EV_(4);
        Z6 = pbm.efpar(23)+pbm.efpar(32)*EV_(5)+pbm.efpar(41)*EV_(5)*EV_(5);
        Z7 = pbm.efpar(24)+pbm.efpar(33)*EV_(7)+pbm.efpar(42)*EV_(7)*EV_(7);
        Z8 = pbm.efpar(25)+pbm.efpar(34)*EV_(8)+pbm.efpar(43)*EV_(8)*EV_(8);
        Z9 = pbm.efpar(26)+pbm.efpar(35)*EV_(9)+pbm.efpar(44)*EV_(9)*EV_(9);
        ZT = pbm.efpar(27)+pbm.efpar(36)*EV_(10)+pbm.efpar(45)*EV_(10)*EV_(10);
        DZ1 = pbm.efpar(28)+2.0*pbm.efpar(37)*EV_(1);
        DZ2 = pbm.efpar(29)+2.0*pbm.efpar(38)*EV_(2);
        DZ3 = pbm.efpar(30)+2.0*pbm.efpar(39)*EV_(3);
        DZ4 = pbm.efpar(31)+2.0*pbm.efpar(40)*EV_(4);
        DZ6 = pbm.efpar(32)+2.0*pbm.efpar(41)*EV_(5);
        DZ7 = pbm.efpar(33)+2.0*pbm.efpar(42)*EV_(7);
        DZ8 = pbm.efpar(34)+2.0*pbm.efpar(43)*EV_(8);
        DZ9 = pbm.efpar(35)+2.0*pbm.efpar(44)*EV_(9);
        DZT = pbm.efpar(36)+2.0*pbm.efpar(45)*EV_(10);
        E1 = pbm.efpar(1)+pbm.efpar(10)*Z1;
        E2 = pbm.efpar(2)+pbm.efpar(11)*(Z2-Z3);
        E3 = pbm.efpar(3)+pbm.efpar(12)*(Z3-Z9);
        E4 = pbm.efpar(4)+pbm.efpar(13)*(Z4-Z9);
        E6 = pbm.efpar(5)+pbm.efpar(14)*(Z6-Z7);
        E7 = pbm.efpar(6)+pbm.efpar(15)*(Z7-Z8);
        E8 = pbm.efpar(7)+pbm.efpar(16)*Z8;
        E9 = pbm.efpar(8)+pbm.efpar(17)*Z9;
        ET = pbm.efpar(9)+pbm.efpar(18)*ZT;
        DE11 = pbm.efpar(10)*DZ1;
        DE22 = pbm.efpar(11)*DZ2;
        DE23 = -pbm.efpar(11)*DZ3;
        DE33 = pbm.efpar(12)*DZ3;
        DE39 = -pbm.efpar(12)*DZ9;
        DE44 = pbm.efpar(13)*DZ4;
        DE49 = -pbm.efpar(13)*DZ9;
        DE66 = pbm.efpar(14)*DZ6;
        DE67 = -pbm.efpar(14)*DZ7;
        DE77 = pbm.efpar(15)*DZ7;
        DE78 = -pbm.efpar(15)*DZ8;
        DE88 = pbm.efpar(16)*DZ8;
        DE99 = pbm.efpar(17)*DZ9;
        DETT = pbm.efpar(18)*DZT;
        HE111 = 2.0*pbm.efpar(10)*pbm.efpar(37);
        HE222 = 2.0*pbm.efpar(11)*pbm.efpar(38);
        HE233 = -2.0*pbm.efpar(11)*pbm.efpar(39);
        HE333 = 2.0*pbm.efpar(12)*pbm.efpar(39);
        HE399 = -2.0*pbm.efpar(12)*pbm.efpar(44);
        HE444 = 2.0*pbm.efpar(13)*pbm.efpar(40);
        HE499 = -2.0*pbm.efpar(13)*pbm.efpar(44);
        HE666 = 2.0*pbm.efpar(14)*pbm.efpar(41);
        HE677 = -2.0*pbm.efpar(14)*pbm.efpar(42);
        HE777 = 2.0*pbm.efpar(15)*pbm.efpar(42);
        HE788 = -2.0*pbm.efpar(15)*pbm.efpar(43);
        HE888 = 2.0*pbm.efpar(16)*pbm.efpar(43);
        HE999 = 2.0*pbm.efpar(17)*pbm.efpar(44);
        HETTT = 2.0*pbm.efpar(18)*pbm.efpar(45);
        PS = EV_(11)*E1+EV_(12)*E2+EV_(13)*E3+EV_(14)*E4+EV_(6)*E6+EV_(15)*...
             E7+EV_(16)*E8+EV_(17)*E9+EV_(18)*ET;
        DP1 = EV_(11)*DE11;
        DP2 = EV_(12)*DE22;
        DP3 = EV_(13)*DE33+EV_(12)*DE23;
        DP4 = EV_(14)*DE44;
        DP6 = EV_(6)*DE66;
        DP7 = EV_(15)*DE77+EV_(6)*DE67;
        DP8 = EV_(16)*DE88+EV_(15)*DE78;
        DP9 = EV_(17)*DE99+EV_(13)*DE39+EV_(14)*DE49;
        DPT = EV_(18)*DETT;
        HP11 = EV_(11)*HE111;
        HP22 = EV_(12)*HE222;
        HP33 = EV_(13)*HE333+EV_(12)*HE233;
        HP44 = EV_(14)*HE444;
        HP66 = EV_(6)*HE666;
        HP77 = EV_(15)*HE777+EV_(6)*HE677;
        HP88 = EV_(16)*HE888+EV_(15)*HE788;
        HP99 = EV_(17)*HE999+EV_(13)*HE399+EV_(14)*HE499;
        HPTT = EV_(18)*HETTT;
        HFPS = pbm.elpar{iel_}(3)+pbm.elpar{iel_}(3);
        DFPS = pbm.elpar{iel_}(2)+HFPS*PS;
        varargout{1} = pbm.elpar{iel_}(1)+pbm.elpar{iel_}(2)*PS+pbm.elpar{iel_}(3)*...
             PS*PS;
        if(nargout>1)
            g_(1,1) = DP1*DFPS;
            g_(2,1) = DP2*DFPS;
            g_(3,1) = DP3*DFPS;
            g_(4,1) = DP4*DFPS;
            g_(5,1) = DP6*DFPS;
            g_(7,1) = DP7*DFPS;
            g_(8,1) = DP8*DFPS;
            g_(9,1) = DP9*DFPS;
            g_(10,1) = DPT*DFPS;
            g_(11,1) = E1*DFPS;
            g_(12,1) = E2*DFPS;
            g_(13,1) = E3*DFPS;
            g_(14,1) = E4*DFPS;
            g_(6,1) = E6*DFPS;
            g_(15,1) = E7*DFPS;
            g_(16,1) = E8*DFPS;
            g_(17,1) = E9*DFPS;
            g_(18,1) = ET*DFPS;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(18,18);
                H_(1,1) = HFPS*DP1*DP1+DFPS*HP11;
                H_(1,2) = HFPS*DP1*DP2;
                H_(2,1) = H_(1,2);
                H_(1,3) = HFPS*DP1*DP3;
                H_(3,1) = H_(1,3);
                H_(1,4) = HFPS*DP1*DP4;
                H_(4,1) = H_(1,4);
                H_(1,5) = HFPS*DP1*DP6;
                H_(5,1) = H_(1,5);
                H_(1,7) = HFPS*DP1*DP7;
                H_(7,1) = H_(1,7);
                H_(1,8) = HFPS*DP1*DP8;
                H_(8,1) = H_(1,8);
                H_(1,9) = HFPS*DP1*DP9;
                H_(9,1) = H_(1,9);
                H_(1,10) = HFPS*DP1*DPT;
                H_(10,1) = H_(1,10);
                H_(1,11) = HFPS*DP1*E1+DFPS*DE11;
                H_(11,1) = H_(1,11);
                H_(1,12) = HFPS*DP1*E2;
                H_(12,1) = H_(1,12);
                H_(1,13) = HFPS*DP1*E3;
                H_(13,1) = H_(1,13);
                H_(1,14) = HFPS*DP1*E4;
                H_(14,1) = H_(1,14);
                H_(1,6) = HFPS*DP1*E6;
                H_(6,1) = H_(1,6);
                H_(1,15) = HFPS*DP1*E7;
                H_(15,1) = H_(1,15);
                H_(1,16) = HFPS*DP1*E8;
                H_(16,1) = H_(1,16);
                H_(1,17) = HFPS*DP1*E9;
                H_(17,1) = H_(1,17);
                H_(1,18) = HFPS*DP1*ET;
                H_(18,1) = H_(1,18);
                H_(2,2) = HFPS*DP2*DP2+DFPS*HP22;
                H_(2,3) = HFPS*DP2*DP3;
                H_(3,2) = H_(2,3);
                H_(2,4) = HFPS*DP2*DP4;
                H_(4,2) = H_(2,4);
                H_(2,5) = HFPS*DP2*DP6;
                H_(5,2) = H_(2,5);
                H_(2,7) = HFPS*DP2*DP7;
                H_(7,2) = H_(2,7);
                H_(2,8) = HFPS*DP2*DP8;
                H_(8,2) = H_(2,8);
                H_(2,9) = HFPS*DP2*DP9;
                H_(9,2) = H_(2,9);
                H_(2,10) = HFPS*DP2*DPT;
                H_(10,2) = H_(2,10);
                H_(2,11) = HFPS*DP2*E1;
                H_(11,2) = H_(2,11);
                H_(2,12) = HFPS*DP2*E2+DFPS*DE22;
                H_(12,2) = H_(2,12);
                H_(2,13) = HFPS*DP2*E3;
                H_(13,2) = H_(2,13);
                H_(2,14) = HFPS*DP2*E4;
                H_(14,2) = H_(2,14);
                H_(2,6) = HFPS*DP2*E6;
                H_(6,2) = H_(2,6);
                H_(2,15) = HFPS*DP2*E7;
                H_(15,2) = H_(2,15);
                H_(2,16) = HFPS*DP2*E8;
                H_(16,2) = H_(2,16);
                H_(2,17) = HFPS*DP2*E9;
                H_(17,2) = H_(2,17);
                H_(2,18) = HFPS*DP2*ET;
                H_(18,2) = H_(2,18);
                H_(3,3) = HFPS*DP3*DP3+DFPS*HP33;
                H_(3,4) = HFPS*DP3*DP4;
                H_(4,3) = H_(3,4);
                H_(3,5) = HFPS*DP3*DP6;
                H_(5,3) = H_(3,5);
                H_(3,7) = HFPS*DP3*DP7;
                H_(7,3) = H_(3,7);
                H_(3,8) = HFPS*DP3*DP8;
                H_(8,3) = H_(3,8);
                H_(3,9) = HFPS*DP3*DP9;
                H_(9,3) = H_(3,9);
                H_(3,10) = HFPS*DP3*DPT;
                H_(10,3) = H_(3,10);
                H_(3,11) = HFPS*DP3*E1;
                H_(11,3) = H_(3,11);
                H_(3,12) = HFPS*DP3*E2+DFPS*DE23;
                H_(12,3) = H_(3,12);
                H_(3,13) = HFPS*DP3*E3+DFPS*DE33;
                H_(13,3) = H_(3,13);
                H_(3,14) = HFPS*DP3*E4;
                H_(14,3) = H_(3,14);
                H_(3,6) = HFPS*DP3*E6;
                H_(6,3) = H_(3,6);
                H_(3,15) = HFPS*DP3*E7;
                H_(15,3) = H_(3,15);
                H_(3,16) = HFPS*DP3*E8;
                H_(16,3) = H_(3,16);
                H_(3,17) = HFPS*DP3*E9;
                H_(17,3) = H_(3,17);
                H_(3,18) = HFPS*DP3*ET;
                H_(18,3) = H_(3,18);
                H_(4,4) = HFPS*DP4*DP4+DFPS*HP44;
                H_(4,5) = HFPS*DP4*DP6;
                H_(5,4) = H_(4,5);
                H_(4,7) = HFPS*DP4*DP7;
                H_(7,4) = H_(4,7);
                H_(4,8) = HFPS*DP4*DP8;
                H_(8,4) = H_(4,8);
                H_(4,9) = HFPS*DP4*DP9;
                H_(9,4) = H_(4,9);
                H_(4,10) = HFPS*DP4*DPT;
                H_(10,4) = H_(4,10);
                H_(4,11) = HFPS*DP4*E1;
                H_(11,4) = H_(4,11);
                H_(4,12) = HFPS*DP4*E2;
                H_(12,4) = H_(4,12);
                H_(4,13) = HFPS*DP4*E3;
                H_(13,4) = H_(4,13);
                H_(4,14) = HFPS*DP4*E4+DFPS*DE44;
                H_(14,4) = H_(4,14);
                H_(4,6) = HFPS*DP4*E6;
                H_(6,4) = H_(4,6);
                H_(4,15) = HFPS*DP4*E7;
                H_(15,4) = H_(4,15);
                H_(4,16) = HFPS*DP4*E8;
                H_(16,4) = H_(4,16);
                H_(4,17) = HFPS*DP4*E9;
                H_(17,4) = H_(4,17);
                H_(4,18) = HFPS*DP4*ET;
                H_(18,4) = H_(4,18);
                H_(5,5) = HFPS*DP6*DP6+DFPS*HP66;
                H_(5,7) = HFPS*DP6*DP7;
                H_(7,5) = H_(5,7);
                H_(5,8) = HFPS*DP6*DP8;
                H_(8,5) = H_(5,8);
                H_(5,9) = HFPS*DP6*DP9;
                H_(9,5) = H_(5,9);
                H_(5,10) = HFPS*DP6*DPT;
                H_(10,5) = H_(5,10);
                H_(5,11) = HFPS*DP6*E1;
                H_(11,5) = H_(5,11);
                H_(5,12) = HFPS*DP6*E2;
                H_(12,5) = H_(5,12);
                H_(5,13) = HFPS*DP6*E3;
                H_(13,5) = H_(5,13);
                H_(5,14) = HFPS*DP6*E4;
                H_(14,5) = H_(5,14);
                H_(5,6) = HFPS*DP6*E6+DFPS*DE66;
                H_(6,5) = H_(5,6);
                H_(5,15) = HFPS*DP6*E7;
                H_(15,5) = H_(5,15);
                H_(5,16) = HFPS*DP6*E8;
                H_(16,5) = H_(5,16);
                H_(5,17) = HFPS*DP6*E9;
                H_(17,5) = H_(5,17);
                H_(5,18) = HFPS*DP6*ET;
                H_(18,5) = H_(5,18);
                H_(7,7) = HFPS*DP7*DP7+DFPS*HP77;
                H_(7,8) = HFPS*DP7*DP8;
                H_(8,7) = H_(7,8);
                H_(7,9) = HFPS*DP7*DP9;
                H_(9,7) = H_(7,9);
                H_(7,10) = HFPS*DP7*DPT;
                H_(10,7) = H_(7,10);
                H_(7,11) = HFPS*DP7*E1;
                H_(11,7) = H_(7,11);
                H_(7,12) = HFPS*DP7*E2;
                H_(12,7) = H_(7,12);
                H_(7,13) = HFPS*DP7*E3;
                H_(13,7) = H_(7,13);
                H_(7,14) = HFPS*DP7*E4;
                H_(14,7) = H_(7,14);
                H_(7,6) = HFPS*DP7*E6+DFPS*DE67;
                H_(6,7) = H_(7,6);
                H_(7,15) = HFPS*DP7*E7+DFPS*DE77;
                H_(15,7) = H_(7,15);
                H_(7,16) = HFPS*DP7*E8;
                H_(16,7) = H_(7,16);
                H_(7,17) = HFPS*DP7*E9;
                H_(17,7) = H_(7,17);
                H_(7,18) = HFPS*DP7*ET;
                H_(18,7) = H_(7,18);
                H_(8,8) = HFPS*DP8*DP8+DFPS*HP88;
                H_(8,9) = HFPS*DP8*DP9;
                H_(9,8) = H_(8,9);
                H_(8,10) = HFPS*DP8*DPT;
                H_(10,8) = H_(8,10);
                H_(8,11) = HFPS*DP8*E1;
                H_(11,8) = H_(8,11);
                H_(8,12) = HFPS*DP8*E2;
                H_(12,8) = H_(8,12);
                H_(8,13) = HFPS*DP8*E3;
                H_(13,8) = H_(8,13);
                H_(8,14) = HFPS*DP8*E4;
                H_(14,8) = H_(8,14);
                H_(8,6) = HFPS*DP8*E6;
                H_(6,8) = H_(8,6);
                H_(8,15) = HFPS*DP8*E7+DFPS*DE78;
                H_(15,8) = H_(8,15);
                H_(8,16) = HFPS*DP8*E8+DFPS*DE88;
                H_(16,8) = H_(8,16);
                H_(8,17) = HFPS*DP8*E9;
                H_(17,8) = H_(8,17);
                H_(8,18) = HFPS*DP8*ET;
                H_(18,8) = H_(8,18);
                H_(9,9) = HFPS*DP9*DP9+DFPS*HP99;
                H_(9,10) = HFPS*DP9*DPT;
                H_(10,9) = H_(9,10);
                H_(9,11) = HFPS*DP9*E1;
                H_(11,9) = H_(9,11);
                H_(9,12) = HFPS*DP9*E2;
                H_(12,9) = H_(9,12);
                H_(9,13) = HFPS*DP9*E3+DFPS*DE39;
                H_(13,9) = H_(9,13);
                H_(9,14) = HFPS*DP9*E4+DFPS*DE49;
                H_(14,9) = H_(9,14);
                H_(9,6) = HFPS*DP9*E6;
                H_(6,9) = H_(9,6);
                H_(9,15) = HFPS*DP9*E7;
                H_(15,9) = H_(9,15);
                H_(9,16) = HFPS*DP9*E8;
                H_(16,9) = H_(9,16);
                H_(9,17) = HFPS*DP9*E9+DFPS*DE99;
                H_(17,9) = H_(9,17);
                H_(9,18) = HFPS*DP9*ET;
                H_(18,9) = H_(9,18);
                H_(10,10) = HFPS*DPT*DPT+DFPS*HPTT;
                H_(10,11) = HFPS*DPT*E1;
                H_(11,10) = H_(10,11);
                H_(10,12) = HFPS*DPT*E2;
                H_(12,10) = H_(10,12);
                H_(10,13) = HFPS*DPT*E3;
                H_(13,10) = H_(10,13);
                H_(10,14) = HFPS*DPT*E4;
                H_(14,10) = H_(10,14);
                H_(10,6) = HFPS*DPT*E6;
                H_(6,10) = H_(10,6);
                H_(10,15) = HFPS*DPT*E7;
                H_(15,10) = H_(10,15);
                H_(10,16) = HFPS*DPT*E8;
                H_(16,10) = H_(10,16);
                H_(10,17) = HFPS*DPT*E9;
                H_(17,10) = H_(10,17);
                H_(10,18) = HFPS*DPT*ET+DFPS*DETT;
                H_(18,10) = H_(10,18);
                H_(11,11) = HFPS*E1*E1;
                H_(11,12) = HFPS*E1*E2;
                H_(12,11) = H_(11,12);
                H_(11,13) = HFPS*E1*E3;
                H_(13,11) = H_(11,13);
                H_(11,14) = HFPS*E1*E4;
                H_(14,11) = H_(11,14);
                H_(11,6) = HFPS*E1*E6;
                H_(6,11) = H_(11,6);
                H_(11,15) = HFPS*E1*E7;
                H_(15,11) = H_(11,15);
                H_(11,16) = HFPS*E1*E8;
                H_(16,11) = H_(11,16);
                H_(11,17) = HFPS*E1*E9;
                H_(17,11) = H_(11,17);
                H_(11,18) = HFPS*E1*ET;
                H_(18,11) = H_(11,18);
                H_(12,12) = HFPS*E2*E2;
                H_(12,13) = HFPS*E2*E3;
                H_(13,12) = H_(12,13);
                H_(12,14) = HFPS*E2*E4;
                H_(14,12) = H_(12,14);
                H_(12,6) = HFPS*E2*E6;
                H_(6,12) = H_(12,6);
                H_(12,15) = HFPS*E2*E7;
                H_(15,12) = H_(12,15);
                H_(12,16) = HFPS*E2*E8;
                H_(16,12) = H_(12,16);
                H_(12,17) = HFPS*E2*E9;
                H_(17,12) = H_(12,17);
                H_(12,18) = HFPS*E2*ET;
                H_(18,12) = H_(12,18);
                H_(13,13) = HFPS*E3*E3;
                H_(13,14) = HFPS*E3*E4;
                H_(14,13) = H_(13,14);
                H_(13,6) = HFPS*E3*E6;
                H_(6,13) = H_(13,6);
                H_(13,15) = HFPS*E3*E7;
                H_(15,13) = H_(13,15);
                H_(13,16) = HFPS*E3*E8;
                H_(16,13) = H_(13,16);
                H_(13,17) = HFPS*E3*E9;
                H_(17,13) = H_(13,17);
                H_(13,18) = HFPS*E3*ET;
                H_(18,13) = H_(13,18);
                H_(14,14) = HFPS*E4*E4;
                H_(14,6) = HFPS*E4*E6;
                H_(6,14) = H_(14,6);
                H_(14,15) = HFPS*E4*E7;
                H_(15,14) = H_(14,15);
                H_(14,16) = HFPS*E4*E8;
                H_(16,14) = H_(14,16);
                H_(14,17) = HFPS*E4*E9;
                H_(17,14) = H_(14,17);
                H_(14,18) = HFPS*E4*ET;
                H_(18,14) = H_(14,18);
                H_(6,6) = HFPS*E6*E6;
                H_(6,15) = HFPS*E6*E7;
                H_(15,6) = H_(6,15);
                H_(6,16) = HFPS*E6*E8;
                H_(16,6) = H_(6,16);
                H_(6,17) = HFPS*E6*E9;
                H_(17,6) = H_(6,17);
                H_(6,18) = HFPS*E6*ET;
                H_(18,6) = H_(6,18);
                H_(15,15) = HFPS*E7*E7;
                H_(15,16) = HFPS*E7*E8;
                H_(16,15) = H_(15,16);
                H_(15,17) = HFPS*E7*E9;
                H_(17,15) = H_(15,17);
                H_(15,18) = HFPS*E7*ET;
                H_(18,15) = H_(15,18);
                H_(16,16) = HFPS*E8*E8;
                H_(16,17) = HFPS*E8*E9;
                H_(17,16) = H_(16,17);
                H_(16,18) = HFPS*E8*ET;
                H_(18,16) = H_(16,18);
                H_(17,17) = HFPS*E9*E9;
                H_(17,18) = HFPS*E9*ET;
                H_(18,17) = H_(17,18);
                H_(18,18) = HFPS*ET*ET;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [45,0];
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

