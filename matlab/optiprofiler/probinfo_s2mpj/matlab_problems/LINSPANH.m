function varargout = LINSPANH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LINSPANH
%    *********
% 
%    A linear network problem based on the spanish hydro-electric
%    reservoir management problem SPANHYD
% 
%    Source:
%    A partial specification of problem SPANHYD.
% 
%    SIF input: Ph. Toint, Sept 1990.
% 
%    classification = 'C-CLNR2-MN-97-33'
% 
%    Number of arcs = 97
%    Number of nodes
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LINSPANH';

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
        ig = ig_('OBJ');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
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
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                77.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'C-CLNR2-MN-97-33';
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

