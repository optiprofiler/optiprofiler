function varargout = SWOPF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    An optimal electrical powerflow system design problem from Switzerland.
% 
%    Source:
%    a contribution to fullfill the LANCELOT academic licence agreement.
% 
%    SIF input: R. Bacher, Dept of Electrical Engineering, ETH Zurich, 
%               November 1994.
% 
%    classification = 'LQR2-RN-83-92'
% 
%    Number of nodes       =   7
%    Number of branches    =   7
%    Number of generators  =   3
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SWOPF';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('FIRST') = 1;
        v_('NOBRANCHES') = 7;
        v_('NOSHUNTS') = 0;
        v_('NOTRAFOS') = 3;
        v_('NOBUSSES') = 7;
        v_('NOGEN') = 3;
        v_('NOGENBK') = 3;
        v_('NOAREAS') = 0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','VE0001',ix_);
        pb.xnames{iv} = 'VE0001';
        [iv,ix_] = s2mpjlib('ii','VF0001',ix_);
        pb.xnames{iv} = 'VF0001';
        [iv,ix_] = s2mpjlib('ii','V20001',ix_);
        pb.xnames{iv} = 'V20001';
        [iv,ix_] = s2mpjlib('ii','VE0002',ix_);
        pb.xnames{iv} = 'VE0002';
        [iv,ix_] = s2mpjlib('ii','VF0002',ix_);
        pb.xnames{iv} = 'VF0002';
        [iv,ix_] = s2mpjlib('ii','V20002',ix_);
        pb.xnames{iv} = 'V20002';
        [iv,ix_] = s2mpjlib('ii','VE0003',ix_);
        pb.xnames{iv} = 'VE0003';
        [iv,ix_] = s2mpjlib('ii','VF0003',ix_);
        pb.xnames{iv} = 'VF0003';
        [iv,ix_] = s2mpjlib('ii','V20003',ix_);
        pb.xnames{iv} = 'V20003';
        [iv,ix_] = s2mpjlib('ii','VE0004',ix_);
        pb.xnames{iv} = 'VE0004';
        [iv,ix_] = s2mpjlib('ii','VF0004',ix_);
        pb.xnames{iv} = 'VF0004';
        [iv,ix_] = s2mpjlib('ii','V20004',ix_);
        pb.xnames{iv} = 'V20004';
        [iv,ix_] = s2mpjlib('ii','VE0005',ix_);
        pb.xnames{iv} = 'VE0005';
        [iv,ix_] = s2mpjlib('ii','VF0005',ix_);
        pb.xnames{iv} = 'VF0005';
        [iv,ix_] = s2mpjlib('ii','V20005',ix_);
        pb.xnames{iv} = 'V20005';
        [iv,ix_] = s2mpjlib('ii','VE0006',ix_);
        pb.xnames{iv} = 'VE0006';
        [iv,ix_] = s2mpjlib('ii','VF0006',ix_);
        pb.xnames{iv} = 'VF0006';
        [iv,ix_] = s2mpjlib('ii','V20006',ix_);
        pb.xnames{iv} = 'V20006';
        [iv,ix_] = s2mpjlib('ii','VE0007',ix_);
        pb.xnames{iv} = 'VE0007';
        [iv,ix_] = s2mpjlib('ii','VF0007',ix_);
        pb.xnames{iv} = 'VF0007';
        [iv,ix_] = s2mpjlib('ii','V20007',ix_);
        pb.xnames{iv} = 'V20007';
        [iv,ix_] = s2mpjlib('ii','EI0001',ix_);
        pb.xnames{iv} = 'EI0001';
        [iv,ix_] = s2mpjlib('ii','FI0001',ix_);
        pb.xnames{iv} = 'FI0001';
        [iv,ix_] = s2mpjlib('ii','EJ0001',ix_);
        pb.xnames{iv} = 'EJ0001';
        [iv,ix_] = s2mpjlib('ii','FJ0001',ix_);
        pb.xnames{iv} = 'FJ0001';
        [iv,ix_] = s2mpjlib('ii','PI0001',ix_);
        pb.xnames{iv} = 'PI0001';
        [iv,ix_] = s2mpjlib('ii','QI0001',ix_);
        pb.xnames{iv} = 'QI0001';
        [iv,ix_] = s2mpjlib('ii','PJ0001',ix_);
        pb.xnames{iv} = 'PJ0001';
        [iv,ix_] = s2mpjlib('ii','QJ0001',ix_);
        pb.xnames{iv} = 'QJ0001';
        [iv,ix_] = s2mpjlib('ii','EI0002',ix_);
        pb.xnames{iv} = 'EI0002';
        [iv,ix_] = s2mpjlib('ii','FI0002',ix_);
        pb.xnames{iv} = 'FI0002';
        [iv,ix_] = s2mpjlib('ii','EJ0002',ix_);
        pb.xnames{iv} = 'EJ0002';
        [iv,ix_] = s2mpjlib('ii','FJ0002',ix_);
        pb.xnames{iv} = 'FJ0002';
        [iv,ix_] = s2mpjlib('ii','PI0002',ix_);
        pb.xnames{iv} = 'PI0002';
        [iv,ix_] = s2mpjlib('ii','QI0002',ix_);
        pb.xnames{iv} = 'QI0002';
        [iv,ix_] = s2mpjlib('ii','PJ0002',ix_);
        pb.xnames{iv} = 'PJ0002';
        [iv,ix_] = s2mpjlib('ii','QJ0002',ix_);
        pb.xnames{iv} = 'QJ0002';
        [iv,ix_] = s2mpjlib('ii','EI0003',ix_);
        pb.xnames{iv} = 'EI0003';
        [iv,ix_] = s2mpjlib('ii','FI0003',ix_);
        pb.xnames{iv} = 'FI0003';
        [iv,ix_] = s2mpjlib('ii','EJ0003',ix_);
        pb.xnames{iv} = 'EJ0003';
        [iv,ix_] = s2mpjlib('ii','FJ0003',ix_);
        pb.xnames{iv} = 'FJ0003';
        [iv,ix_] = s2mpjlib('ii','PI0003',ix_);
        pb.xnames{iv} = 'PI0003';
        [iv,ix_] = s2mpjlib('ii','QI0003',ix_);
        pb.xnames{iv} = 'QI0003';
        [iv,ix_] = s2mpjlib('ii','PJ0003',ix_);
        pb.xnames{iv} = 'PJ0003';
        [iv,ix_] = s2mpjlib('ii','QJ0003',ix_);
        pb.xnames{iv} = 'QJ0003';
        [iv,ix_] = s2mpjlib('ii','EI0004',ix_);
        pb.xnames{iv} = 'EI0004';
        [iv,ix_] = s2mpjlib('ii','FI0004',ix_);
        pb.xnames{iv} = 'FI0004';
        [iv,ix_] = s2mpjlib('ii','EJ0004',ix_);
        pb.xnames{iv} = 'EJ0004';
        [iv,ix_] = s2mpjlib('ii','FJ0004',ix_);
        pb.xnames{iv} = 'FJ0004';
        [iv,ix_] = s2mpjlib('ii','PI0004',ix_);
        pb.xnames{iv} = 'PI0004';
        [iv,ix_] = s2mpjlib('ii','QI0004',ix_);
        pb.xnames{iv} = 'QI0004';
        [iv,ix_] = s2mpjlib('ii','PJ0004',ix_);
        pb.xnames{iv} = 'PJ0004';
        [iv,ix_] = s2mpjlib('ii','QJ0004',ix_);
        pb.xnames{iv} = 'QJ0004';
        [iv,ix_] = s2mpjlib('ii','EI0005',ix_);
        pb.xnames{iv} = 'EI0005';
        [iv,ix_] = s2mpjlib('ii','FI0005',ix_);
        pb.xnames{iv} = 'FI0005';
        [iv,ix_] = s2mpjlib('ii','EJ0005',ix_);
        pb.xnames{iv} = 'EJ0005';
        [iv,ix_] = s2mpjlib('ii','FJ0005',ix_);
        pb.xnames{iv} = 'FJ0005';
        [iv,ix_] = s2mpjlib('ii','PI0005',ix_);
        pb.xnames{iv} = 'PI0005';
        [iv,ix_] = s2mpjlib('ii','QI0005',ix_);
        pb.xnames{iv} = 'QI0005';
        [iv,ix_] = s2mpjlib('ii','PJ0005',ix_);
        pb.xnames{iv} = 'PJ0005';
        [iv,ix_] = s2mpjlib('ii','QJ0005',ix_);
        pb.xnames{iv} = 'QJ0005';
        [iv,ix_] = s2mpjlib('ii','EI0006',ix_);
        pb.xnames{iv} = 'EI0006';
        [iv,ix_] = s2mpjlib('ii','FI0006',ix_);
        pb.xnames{iv} = 'FI0006';
        [iv,ix_] = s2mpjlib('ii','EJ0006',ix_);
        pb.xnames{iv} = 'EJ0006';
        [iv,ix_] = s2mpjlib('ii','FJ0006',ix_);
        pb.xnames{iv} = 'FJ0006';
        [iv,ix_] = s2mpjlib('ii','PI0006',ix_);
        pb.xnames{iv} = 'PI0006';
        [iv,ix_] = s2mpjlib('ii','QI0006',ix_);
        pb.xnames{iv} = 'QI0006';
        [iv,ix_] = s2mpjlib('ii','PJ0006',ix_);
        pb.xnames{iv} = 'PJ0006';
        [iv,ix_] = s2mpjlib('ii','QJ0006',ix_);
        pb.xnames{iv} = 'QJ0006';
        [iv,ix_] = s2mpjlib('ii','EI0007',ix_);
        pb.xnames{iv} = 'EI0007';
        [iv,ix_] = s2mpjlib('ii','FI0007',ix_);
        pb.xnames{iv} = 'FI0007';
        [iv,ix_] = s2mpjlib('ii','EJ0007',ix_);
        pb.xnames{iv} = 'EJ0007';
        [iv,ix_] = s2mpjlib('ii','FJ0007',ix_);
        pb.xnames{iv} = 'FJ0007';
        [iv,ix_] = s2mpjlib('ii','PI0007',ix_);
        pb.xnames{iv} = 'PI0007';
        [iv,ix_] = s2mpjlib('ii','QI0007',ix_);
        pb.xnames{iv} = 'QI0007';
        [iv,ix_] = s2mpjlib('ii','PJ0007',ix_);
        pb.xnames{iv} = 'PJ0007';
        [iv,ix_] = s2mpjlib('ii','QJ0007',ix_);
        pb.xnames{iv} = 'QJ0007';
        [iv,ix_] = s2mpjlib('ii','PG0001',ix_);
        pb.xnames{iv} = 'PG0001';
        [iv,ix_] = s2mpjlib('ii','PG0002',ix_);
        pb.xnames{iv} = 'PG0002';
        [iv,ix_] = s2mpjlib('ii','PG0003',ix_);
        pb.xnames{iv} = 'PG0003';
        [iv,ix_] = s2mpjlib('ii','QG0001',ix_);
        pb.xnames{iv} = 'QG0001';
        [iv,ix_] = s2mpjlib('ii','QG0002',ix_);
        pb.xnames{iv} = 'QG0002';
        [iv,ix_] = s2mpjlib('ii','QG0003',ix_);
        pb.xnames{iv} = 'QG0003';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','GV20001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20001';
        iv = ix_('V20001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','SLF0000',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SLF0000';
        iv = ix_('VF0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20002';
        iv = ix_('V20002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20003';
        iv = ix_('V20003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20004';
        iv = ix_('V20004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20005';
        iv = ix_('V20005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20006';
        iv = ix_('V20006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GV20007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GV20007';
        iv = ix_('V20007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0001';
        iv = ix_('EI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0001';
        iv = ix_('FI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0001';
        iv = ix_('EJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0001';
        iv = ix_('FJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0001';
        iv = ix_('VE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.299;
        end
        iv = ix_('VF0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -66.243;
        end
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.299;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.243;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0001';
        iv = ix_('VE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.243;
        end
        iv = ix_('VF0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.299;
        end
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -66.243;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.299;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0001';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.299;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -66.243;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0001';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.243;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.299;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0001';
        iv = ix_('VE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.299;
        end
        iv = ix_('VF0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.243;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0001';
        iv = ix_('VE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -66.243+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -66.243;
        end
        iv = ix_('VF0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.299+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.299;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0001';
        iv = ix_('PI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0001';
        iv = ix_('QI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0001';
        iv = ix_('PJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0001';
        iv = ix_('QJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0001';
        iv = ix_('PI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0002';
        iv = ix_('PJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0001';
        iv = ix_('QI0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0002';
        iv = ix_('QJ0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0001';
        [ig,ig_] = s2mpjlib('ii','GMXJ0001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0001';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0002';
        iv = ix_('EI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0002';
        iv = ix_('FI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0002';
        iv = ix_('EJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0002';
        iv = ix_('FJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0002';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.915;
        end
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.051;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0002';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.915;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.051;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0002';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.915;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0002';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.915;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0002';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.051;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0002';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.051;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0002';
        iv = ix_('PI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0002';
        iv = ix_('QI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0002';
        iv = ix_('PJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0002';
        iv = ix_('QJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0002';
        iv = ix_('PI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0006';
        iv = ix_('PJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0002';
        iv = ix_('QI0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0006';
        iv = ix_('QJ0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0002';
        [ig,ig_] = s2mpjlib('ii','GMXJ0002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0002';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0003';
        iv = ix_('EI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0003';
        iv = ix_('FI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0003';
        iv = ix_('EJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0003';
        iv = ix_('FJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0003';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.498;
        end
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.588;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0003';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.498;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.588;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0003';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.498;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0003';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.498;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0003';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.588;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0003';
        iv = ix_('VE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.588;
        end
        iv = ix_('VF0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0003';
        iv = ix_('PI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0003';
        iv = ix_('QI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0003';
        iv = ix_('PJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0003';
        iv = ix_('QJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0002';
        iv = ix_('PI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0004';
        iv = ix_('PJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0002';
        iv = ix_('QI0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0004';
        iv = ix_('QJ0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0003';
        [ig,ig_] = s2mpjlib('ii','GMXJ0003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0003';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0004';
        iv = ix_('EI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0004';
        iv = ix_('FI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0004';
        iv = ix_('EJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0004';
        iv = ix_('FJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0004';
        iv = ix_('VE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.897;
        end
        iv = ix_('VF0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -82.759;
        end
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.897;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 82.759;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0004';
        iv = ix_('VE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 82.759;
        end
        iv = ix_('VF0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.897;
        end
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -82.759;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.897;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0004';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.897;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -82.759;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0004';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 82.759;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.897;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0004';
        iv = ix_('VE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.897;
        end
        iv = ix_('VF0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 82.759;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0004';
        iv = ix_('VE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -82.759+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -82.759;
        end
        iv = ix_('VF0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.897+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.897;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0004';
        iv = ix_('PI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0004';
        iv = ix_('QI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0004';
        iv = ix_('PJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0004';
        iv = ix_('QJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0003';
        iv = ix_('PI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0004';
        iv = ix_('PJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0003';
        iv = ix_('QI0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0004';
        iv = ix_('QJ0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0004';
        [ig,ig_] = s2mpjlib('ii','GMXJ0004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0004';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0005';
        iv = ix_('EI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0005';
        iv = ix_('FI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0005';
        iv = ix_('EJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0005';
        iv = ix_('FJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0005';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.915;
        end
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.051;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0005';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.915;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.051;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0005';
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.915;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0005';
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.915+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.915;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0005';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.051;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0005';
        iv = ix_('VE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.051+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.051;
        end
        iv = ix_('VF0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.175;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0005';
        iv = ix_('PI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0005';
        iv = ix_('QI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0005';
        iv = ix_('PJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0005';
        iv = ix_('QJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0004';
        iv = ix_('PI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0007';
        iv = ix_('PJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0004';
        iv = ix_('QI0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0007';
        iv = ix_('QJ0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0005';
        [ig,ig_] = s2mpjlib('ii','GMXJ0005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0005';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0006';
        iv = ix_('EI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0006';
        iv = ix_('FI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0006';
        iv = ix_('EJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0006';
        iv = ix_('FJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0006';
        iv = ix_('VE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.448;
        end
        iv = ix_('VF0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -41.379;
        end
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.448;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 41.379;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0006';
        iv = ix_('VE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 41.379;
        end
        iv = ix_('VF0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.448;
        end
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -41.379;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.448;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0006';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.448;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -41.379;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0006';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 41.379;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.448;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0006';
        iv = ix_('VE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.448;
        end
        iv = ix_('VF0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 41.379;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0006';
        iv = ix_('VE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -41.379+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -41.379;
        end
        iv = ix_('VF0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.448;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0006';
        iv = ix_('PI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0006';
        iv = ix_('QI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0006';
        iv = ix_('PJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0006';
        iv = ix_('QJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0005';
        iv = ix_('PI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0006';
        iv = ix_('PJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0005';
        iv = ix_('QI0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0006';
        iv = ix_('QJ0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0006';
        [ig,ig_] = s2mpjlib('ii','GMXJ0006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0006';
        [ig,ig_] = s2mpjlib('ii','LOSS0000',ig_);
        gtype{ig} = '<>';
        iv = ix_('PI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        iv = ix_('PJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0007';
        iv = ix_('EI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0007';
        iv = ix_('FI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0007';
        iv = ix_('EJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0007';
        iv = ix_('FJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GEI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEI0007';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.498;
        end
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.588;
        end
        [ig,ig_] = s2mpjlib('ii','GFI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFI0007';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.498;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.588;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0007';
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.498;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0007';
        iv = ix_('VE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.498+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.498;
        end
        iv = ix_('VF0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GEJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GEJ0007';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.588;
        end
        [ig,ig_] = s2mpjlib('ii','GFJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GFJ0007';
        iv = ix_('VE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.588;
        end
        iv = ix_('VF0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.726+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.726;
        end
        [ig,ig_] = s2mpjlib('ii','GPI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPI0007';
        iv = ix_('PI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQI0007';
        iv = ix_('QI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPJ0007';
        iv = ix_('PJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQJ0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQJ0007';
        iv = ix_('QJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0006';
        iv = ix_('PI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0007';
        iv = ix_('PJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0006';
        iv = ix_('QI0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0007';
        iv = ix_('QJ0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GMXI0007',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXI0007';
        [ig,ig_] = s2mpjlib('ii','GMXJ0007',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'GMXJ0007';
        [ig,ig_] = s2mpjlib('ii','GPNI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0001';
        iv = ix_('PG0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0003';
        iv = ix_('PG0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GPNI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GPNI0005';
        iv = ix_('PG0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0001';
        iv = ix_('QG0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0003';
        iv = ix_('QG0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        [ig,ig_] = s2mpjlib('ii','GQNI0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'GQNI0005';
        iv = ix_('QG0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        pbm.gconst(ig_('GPNI0001')) = 0.000;
        pbm.gconst(ig_('GQNI0001')) = 0.000;
        pbm.gconst(ig_('SLF0000')) = 0.000;
        pbm.gconst(ig_('GPNI0002')) = 2.000;
        pbm.gconst(ig_('GQNI0002')) = 3.000;
        pbm.gconst(ig_('GPNI0003')) = 0.600;
        pbm.gconst(ig_('GQNI0003')) = 0.080;
        pbm.gconst(ig_('GPNI0004')) = 2.000;
        pbm.gconst(ig_('GQNI0004')) = 0.200;
        pbm.gconst(ig_('GPNI0005')) = 0.500;
        pbm.gconst(ig_('GQNI0005')) = 0.050;
        pbm.gconst(ig_('GPNI0006')) = 1.000;
        pbm.gconst(ig_('GQNI0006')) = 0.300;
        pbm.gconst(ig_('GPNI0007')) = 2.000;
        pbm.gconst(ig_('GQNI0007')) = 1.000;
        pbm.gconst(ig_('GMXI0001')) = 16.000;
        pbm.gconst(ig_('GMXJ0001')) = 16.000;
        pbm.gconst(ig_('GMXI0002')) = 4.000;
        pbm.gconst(ig_('GMXJ0002')) = 4.000;
        pbm.gconst(ig_('GMXI0003')) = 4.000;
        pbm.gconst(ig_('GMXJ0003')) = 4.000;
        pbm.gconst(ig_('GMXI0004')) = 25.000;
        pbm.gconst(ig_('GMXJ0004')) = 25.000;
        pbm.gconst(ig_('GMXI0005')) = 4.000;
        pbm.gconst(ig_('GMXJ0005')) = 4.000;
        pbm.gconst(ig_('GMXI0006')) = 6.250;
        pbm.gconst(ig_('GMXJ0006')) = 6.250;
        pbm.gconst(ig_('GMXI0007')) = 4.000;
        pbm.gconst(ig_('GMXJ0007')) = 4.000;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('V20001'),1) = 0.810;
        pb.xupper(ix_('V20001')) = 1.210;
        pb.xlower(ix_('V20002'),1) = 0.810;
        pb.xupper(ix_('V20002')) = 1.210;
        pb.xlower(ix_('V20003'),1) = 0.941;
        pb.xupper(ix_('V20003')) = 1.210;
        pb.xlower(ix_('V20004'),1) = 0.941;
        pb.xupper(ix_('V20004')) = 1.210;
        pb.xlower(ix_('V20005'),1) = 0.941;
        pb.xupper(ix_('V20005')) = 1.210;
        pb.xlower(ix_('V20006'),1) = 0.941;
        pb.xupper(ix_('V20006')) = 1.210;
        pb.xlower(ix_('V20007'),1) = 0.941;
        pb.xupper(ix_('V20007')) = 1.210;
        pb.xlower(ix_('PG0001'),1) = 0.500;
        pb.xupper(ix_('PG0001')) = 10.000;
        pb.xlower(ix_('PG0002'),1) = 0.500;
        pb.xupper(ix_('PG0002')) = 10.000;
        pb.xlower(ix_('PG0003'),1) = 0.200;
        pb.xupper(ix_('PG0003')) = 4.000;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('VE0001'),1) = 1.000;
        pb.x0(ix_('VF0001'),1) = 0.000;
        pb.x0(ix_('V20001'),1) = 1.000;
        pb.x0(ix_('VE0002'),1) = 1.001;
        pb.x0(ix_('VF0002'),1) = 0.000;
        pb.x0(ix_('V20002'),1) = 1.002;
        pb.x0(ix_('VE0003'),1) = 1.050;
        pb.x0(ix_('VF0003'),1) = 0.000;
        pb.x0(ix_('V20003'),1) = 1.102;
        pb.x0(ix_('VE0004'),1) = 1.001;
        pb.x0(ix_('VF0004'),1) = 0.000;
        pb.x0(ix_('V20004'),1) = 1.002;
        pb.x0(ix_('VE0005'),1) = 1.050;
        pb.x0(ix_('VF0005'),1) = 0.000;
        pb.x0(ix_('V20005'),1) = 1.102;
        pb.x0(ix_('VE0006'),1) = 1.001;
        pb.x0(ix_('VF0006'),1) = 0.000;
        pb.x0(ix_('V20006'),1) = 1.002;
        pb.x0(ix_('VE0007'),1) = 1.001;
        pb.x0(ix_('VF0007'),1) = 0.000;
        pb.x0(ix_('V20007'),1) = 1.002;
        pb.x0(ix_('EI0001'),1) = -0.005;
        pb.x0(ix_('FI0001'),1) = 0.066;
        pb.x0(ix_('EJ0001'),1) = 0.005;
        pb.x0(ix_('FJ0001'),1) = -0.066;
        pb.x0(ix_('PI0001'),1) = -0.005;
        pb.x0(ix_('QI0001'),1) = -0.066;
        pb.x0(ix_('PJ0001'),1) = 0.005;
        pb.x0(ix_('QJ0001'),1) = 0.066;
        pb.x0(ix_('EI0002'),1) = 0.000;
        pb.x0(ix_('FI0002'),1) = 0.136;
        pb.x0(ix_('EJ0002'),1) = 0.000;
        pb.x0(ix_('FJ0002'),1) = 0.136;
        pb.x0(ix_('PI0002'),1) = 0.000;
        pb.x0(ix_('QI0002'),1) = -0.136;
        pb.x0(ix_('PJ0002'),1) = 0.000;
        pb.x0(ix_('QJ0002'),1) = -0.136;
        pb.x0(ix_('EI0003'),1) = 0.000;
        pb.x0(ix_('FI0003'),1) = 0.091;
        pb.x0(ix_('EJ0003'),1) = 0.000;
        pb.x0(ix_('FJ0003'),1) = 0.091;
        pb.x0(ix_('PI0003'),1) = 0.000;
        pb.x0(ix_('QI0003'),1) = -0.091;
        pb.x0(ix_('PJ0003'),1) = 0.000;
        pb.x0(ix_('QJ0003'),1) = -0.091;
        pb.x0(ix_('EI0004'),1) = 0.338;
        pb.x0(ix_('FI0004'),1) = -4.055;
        pb.x0(ix_('EJ0004'),1) = -0.338;
        pb.x0(ix_('FJ0004'),1) = 4.055;
        pb.x0(ix_('PI0004'),1) = 0.355;
        pb.x0(ix_('QI0004'),1) = 4.258;
        pb.x0(ix_('PJ0004'),1) = -0.338;
        pb.x0(ix_('QJ0004'),1) = -4.059;
        pb.x0(ix_('EI0005'),1) = 0.000;
        pb.x0(ix_('FI0005'),1) = 0.136;
        pb.x0(ix_('EJ0005'),1) = 0.000;
        pb.x0(ix_('FJ0005'),1) = 0.136;
        pb.x0(ix_('PI0005'),1) = 0.000;
        pb.x0(ix_('QI0005'),1) = -0.136;
        pb.x0(ix_('PJ0005'),1) = 0.000;
        pb.x0(ix_('QJ0005'),1) = -0.136;
        pb.x0(ix_('EI0006'),1) = 0.169;
        pb.x0(ix_('FI0006'),1) = -2.028;
        pb.x0(ix_('EJ0006'),1) = -0.169;
        pb.x0(ix_('FJ0006'),1) = 2.028;
        pb.x0(ix_('PI0006'),1) = 0.177;
        pb.x0(ix_('QI0006'),1) = 2.129;
        pb.x0(ix_('PJ0006'),1) = -0.169;
        pb.x0(ix_('QJ0006'),1) = -2.030;
        pb.x0(ix_('EI0007'),1) = 0.000;
        pb.x0(ix_('FI0007'),1) = 0.091;
        pb.x0(ix_('EJ0007'),1) = 0.000;
        pb.x0(ix_('FJ0007'),1) = 0.091;
        pb.x0(ix_('PI0007'),1) = 0.000;
        pb.x0(ix_('QI0007'),1) = -0.091;
        pb.x0(ix_('PJ0007'),1) = 0.000;
        pb.x0(ix_('QJ0007'),1) = -0.091;
        pb.x0(ix_('PG0001'),1) = 3.000;
        pb.x0(ix_('PG0002'),1) = 5.000;
        pb.x0(ix_('PG0003'),1) = 2.000;
        pb.x0(ix_('QG0001'),1) = 0.000;
        pb.x0(ix_('QG0002'),1) = 0.000;
        pb.x0(ix_('QG0003'),1) = 0.000;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXTIMESY',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eXSQUARE',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VE0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'VF0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20001';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20002';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20003';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20004';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20005';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20006';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIEI0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIFI0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EIFI0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FIEI0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJEJ0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJFJ0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'EJFJ0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'EJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VF0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'FJEJ0007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXTIMESY';
        ielftype(ie) = iet_('eXTIMESY');
        vname = 'FJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'VE0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PI20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QI20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QI0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'PJ20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'PJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'QJ20007';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXSQUARE';
        ielftype(ie) = iet_('eXSQUARE');
        vname = 'QJ0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('GV20001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GV20007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GPI0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20001');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20002');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20003');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20004');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20005');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20006');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPI0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIEI0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIFI0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQI0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EIFI0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FIEI0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GPJ0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJEJ0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJFJ0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        ig = ig_('GQJ0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('EJFJ0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('FJEJ0007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXI0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PI20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QI20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        ig = ig_('GMXJ0007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PJ20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QJ20007');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.000;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solutions
% LO SOLTN               6.78605619D-2
% LO SOLTN               6.78422730D-2
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-RN-83-92';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eXTIMESY'

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

    case 'eXSQUARE'

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
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

