function varargout = CORE1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CORE1
%    *********
% 
%    A problem from the exploitation of a gas transmission network
%    with consideration of the head loses in the pipes. The aim is
%    to satisfy the demand at several points in the network at a
%    minimal pressure, pumping the gas from a number of different
%    entry points.
% 
%    Sources:
%    D. De Wolf, "Optimisation de reseaux de transport de gas avec
%                 consideration des pertes de charge dans les gazoducs",
%                Ph. D. dissertation, CORE, Belgium, 1992, and
%    D. De Wolf, O. Janssens de Bisthoven and Y. Smeers,
%                "The simplex algorithm extended to piecewise linearly
%                 constrained problems II; an application to the gas
%                 transmission problem", CORE discussion paper 9103, 1991.
% 
% 
%    SDIF input: E. Loute and D. De Wolf, September 1992.
% 
%    classification = 'LQI2-RN-65-59'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CORE1';

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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','COST',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','NODE0001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0001';
        [ig,ig_] = s2mpjlib('ii','NODE0002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0002';
        [ig,ig_] = s2mpjlib('ii','NODE0003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0003';
        [ig,ig_] = s2mpjlib('ii','NODE0004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0004';
        [ig,ig_] = s2mpjlib('ii','NODE0005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0005';
        [ig,ig_] = s2mpjlib('ii','NODE0006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0006';
        [ig,ig_] = s2mpjlib('ii','NODE0007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0007';
        [ig,ig_] = s2mpjlib('ii','NODE0008',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0008';
        [ig,ig_] = s2mpjlib('ii','NODE0009',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0009';
        [ig,ig_] = s2mpjlib('ii','NODE0010',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0010';
        [ig,ig_] = s2mpjlib('ii','NODE0011',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0011';
        [ig,ig_] = s2mpjlib('ii','NODE0012',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0012';
        [ig,ig_] = s2mpjlib('ii','NODE0013',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0013';
        [ig,ig_] = s2mpjlib('ii','NODE0014',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0014';
        [ig,ig_] = s2mpjlib('ii','NODE0015',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0015';
        [ig,ig_] = s2mpjlib('ii','NODE0016',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0016';
        [ig,ig_] = s2mpjlib('ii','NODE0017',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0017';
        [ig,ig_] = s2mpjlib('ii','NODE0018',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0018';
        [ig,ig_] = s2mpjlib('ii','NODE0019',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0019';
        [ig,ig_] = s2mpjlib('ii','NODE0020',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'NODE0020';
        [ig,ig_] = s2mpjlib('ii','ARC00001',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00001';
        [ig,ig_] = s2mpjlib('ii','ARC00002',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00002';
        [ig,ig_] = s2mpjlib('ii','ARC00003',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00003';
        [ig,ig_] = s2mpjlib('ii','ARC00004',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00004';
        [ig,ig_] = s2mpjlib('ii','ARC00005',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00005';
        [ig,ig_] = s2mpjlib('ii','ARC00006',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00006';
        [ig,ig_] = s2mpjlib('ii','ARC00007',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00007';
        [ig,ig_] = s2mpjlib('ii','ARC00008',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00008';
        [ig,ig_] = s2mpjlib('ii','ARC00009',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00009';
        [ig,ig_] = s2mpjlib('ii','ARC00010',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'ARC00010';
        [ig,ig_] = s2mpjlib('ii','ARC00011',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'ARC00011';
        [ig,ig_] = s2mpjlib('ii','ARC00012',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00012';
        [ig,ig_] = s2mpjlib('ii','ARC00013',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00013';
        [ig,ig_] = s2mpjlib('ii','ARC00014',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00014';
        [ig,ig_] = s2mpjlib('ii','ARC00015',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00015';
        [ig,ig_] = s2mpjlib('ii','ARC00016',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00016';
        [ig,ig_] = s2mpjlib('ii','ARC00017',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00017';
        [ig,ig_] = s2mpjlib('ii','ARC00018',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00018';
        [ig,ig_] = s2mpjlib('ii','ARC00019',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00019';
        [ig,ig_] = s2mpjlib('ii','ARC00020',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00020';
        [ig,ig_] = s2mpjlib('ii','ARC00021',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00021';
        [ig,ig_] = s2mpjlib('ii','ARC00022',ig_);
        gtype{ig}  = '>=';
        cnames{ig} = 'ARC00022';
        [ig,ig_] = s2mpjlib('ii','ARC00023',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00023';
        [ig,ig_] = s2mpjlib('ii','ARC00024',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'ARC00024';
        [ig,ig_] = s2mpjlib('ii','REGIO001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO001';
        [ig,ig_] = s2mpjlib('ii','REGIO002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO002';
        [ig,ig_] = s2mpjlib('ii','REGIO003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO003';
        [ig,ig_] = s2mpjlib('ii','REGIO004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO004';
        [ig,ig_] = s2mpjlib('ii','REGIO005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO005';
        [ig,ig_] = s2mpjlib('ii','REGIO006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO006';
        [ig,ig_] = s2mpjlib('ii','REGIO007',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO007';
        [ig,ig_] = s2mpjlib('ii','REGIO008',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO008';
        [ig,ig_] = s2mpjlib('ii','REGIO009',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'REGIO009';
        [ig,ig_] = s2mpjlib('ii','PROD0001',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0001';
        [ig,ig_] = s2mpjlib('ii','PROD0002',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0002';
        [ig,ig_] = s2mpjlib('ii','PROD0003',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0003';
        [ig,ig_] = s2mpjlib('ii','PROD0004',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0004';
        [ig,ig_] = s2mpjlib('ii','PROD0005',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0005';
        [ig,ig_] = s2mpjlib('ii','PROD0006',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'PROD0006';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = numEntries(ig_);
        [iv,ix_] = s2mpjlib('ii','FLOW0001',ix_);
        pb.xnames{iv} = 'FLOW0001';
        ig = ig_('NODE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0002',ix_);
        pb.xnames{iv} = 'FLOW0002';
        ig = ig_('NODE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0003',ix_);
        pb.xnames{iv} = 'FLOW0003';
        ig = ig_('NODE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0004',ix_);
        pb.xnames{iv} = 'FLOW0004';
        ig = ig_('NODE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0005',ix_);
        pb.xnames{iv} = 'FLOW0005';
        ig = ig_('NODE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0006',ix_);
        pb.xnames{iv} = 'FLOW0006';
        ig = ig_('NODE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0007',ix_);
        pb.xnames{iv} = 'FLOW0007';
        ig = ig_('NODE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0008',ix_);
        pb.xnames{iv} = 'FLOW0008';
        ig = ig_('NODE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0009',ix_);
        pb.xnames{iv} = 'FLOW0009';
        ig = ig_('NODE0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0010',ix_);
        pb.xnames{iv} = 'FLOW0010';
        ig = ig_('NODE0008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0011',ix_);
        pb.xnames{iv} = 'FLOW0011';
        ig = ig_('NODE0008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0012',ix_);
        pb.xnames{iv} = 'FLOW0012';
        ig = ig_('NODE0009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0013',ix_);
        pb.xnames{iv} = 'FLOW0013';
        ig = ig_('NODE0009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0014',ix_);
        pb.xnames{iv} = 'FLOW0014';
        ig = ig_('NODE0010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0015',ix_);
        pb.xnames{iv} = 'FLOW0015';
        ig = ig_('NODE0010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0016',ix_);
        pb.xnames{iv} = 'FLOW0016';
        ig = ig_('NODE0011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0012');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0017',ix_);
        pb.xnames{iv} = 'FLOW0017';
        ig = ig_('NODE0012');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0013');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0018',ix_);
        pb.xnames{iv} = 'FLOW0018';
        ig = ig_('NODE0013');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0019',ix_);
        pb.xnames{iv} = 'FLOW0019';
        ig = ig_('NODE0014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0015');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0020',ix_);
        pb.xnames{iv} = 'FLOW0020';
        ig = ig_('NODE0015');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0016');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0021',ix_);
        pb.xnames{iv} = 'FLOW0021';
        ig = ig_('NODE0011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0017');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0022',ix_);
        pb.xnames{iv} = 'FLOW0022';
        ig = ig_('NODE0017');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0018');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0023',ix_);
        pb.xnames{iv} = 'FLOW0023';
        ig = ig_('NODE0018');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0019');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','FLOW0024',ix_);
        pb.xnames{iv} = 'FLOW0024';
        ig = ig_('NODE0019');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0020');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00001',ix_);
        pb.xnames{iv} = 'DEM00001';
        ig = ig_('NODE0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00002',ix_);
        pb.xnames{iv} = 'DEM00002';
        ig = ig_('NODE0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00003',ix_);
        pb.xnames{iv} = 'DEM00003';
        ig = ig_('NODE0007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00004',ix_);
        pb.xnames{iv} = 'DEM00004';
        ig = ig_('NODE0010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00005',ix_);
        pb.xnames{iv} = 'DEM00005';
        ig = ig_('NODE0012');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00006',ix_);
        pb.xnames{iv} = 'DEM00006';
        ig = ig_('NODE0015');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00007',ix_);
        pb.xnames{iv} = 'DEM00007';
        ig = ig_('NODE0016');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00008',ix_);
        pb.xnames{iv} = 'DEM00008';
        ig = ig_('NODE0019');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','DEM00009',ix_);
        pb.xnames{iv} = 'DEM00009';
        ig = ig_('NODE0020');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('REGIO009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0001',ix_);
        pb.xnames{iv} = 'SUPP0001';
        ig = ig_('PROD0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0002',ix_);
        pb.xnames{iv} = 'SUPP0002';
        ig = ig_('PROD0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0003',ix_);
        pb.xnames{iv} = 'SUPP0003';
        ig = ig_('PROD0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0004',ix_);
        pb.xnames{iv} = 'SUPP0004';
        ig = ig_('PROD0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0005',ix_);
        pb.xnames{iv} = 'SUPP0005';
        ig = ig_('PROD0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0013');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','SUPP0006',ix_);
        pb.xnames{iv} = 'SUPP0006';
        ig = ig_('PROD0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000E+00;
        end
        ig = ig_('NODE0014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0001',ix_);
        pb.xnames{iv} = 'PROD0001';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.28000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.28000E+00;
        end
        ig = ig_('PROD0001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0002',ix_);
        pb.xnames{iv} = 'PROD0002';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.28000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.28000E+00;
        end
        ig = ig_('PROD0002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0003',ix_);
        pb.xnames{iv} = 'PROD0003';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.28000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.28000E+00;
        end
        ig = ig_('PROD0003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0004',ix_);
        pb.xnames{iv} = 'PROD0004';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.68000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.68000E+00;
        end
        ig = ig_('PROD0004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0005',ix_);
        pb.xnames{iv} = 'PROD0005';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.68000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.68000E+00;
        end
        ig = ig_('PROD0005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PROD0006',ix_);
        pb.xnames{iv} = 'PROD0006';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.68000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.68000E+00;
        end
        ig = ig_('PROD0006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000001',ix_);
        pb.xnames{iv} = 'PI000001';
        ig = ig_('ARC00001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.07027E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.07027E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000001',ix_);
        pb.xnames{iv} = 'PI000001';
        ig = ig_('ARC00002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.07027E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.07027E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000002',ix_);
        pb.xnames{iv} = 'PI000002';
        ig = ig_('ARC00001');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.07027E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.07027E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000002',ix_);
        pb.xnames{iv} = 'PI000002';
        ig = ig_('ARC00002');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.07027E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.07027E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000002',ix_);
        pb.xnames{iv} = 'PI000002';
        ig = ig_('ARC00003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.04685E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.04685E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000002',ix_);
        pb.xnames{iv} = 'PI000002';
        ig = ig_('ARC00004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.04685E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.04685E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000003',ix_);
        pb.xnames{iv} = 'PI000003';
        ig = ig_('ARC00003');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.04685E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.04685E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000003',ix_);
        pb.xnames{iv} = 'PI000003';
        ig = ig_('ARC00004');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.04685E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.04685E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000003',ix_);
        pb.xnames{iv} = 'PI000003';
        ig = ig_('ARC00005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.39543E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.39543E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000004',ix_);
        pb.xnames{iv} = 'PI000004';
        ig = ig_('ARC00005');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.39543E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.39543E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000004',ix_);
        pb.xnames{iv} = 'PI000004';
        ig = ig_('ARC00008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.26895E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.26895E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000004',ix_);
        pb.xnames{iv} = 'PI000004';
        ig = ig_('ARC00009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.59656E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.59656E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000005',ix_);
        pb.xnames{iv} = 'PI000005';
        ig = ig_('ARC00006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00256E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00256E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000006',ix_);
        pb.xnames{iv} = 'PI000006';
        ig = ig_('ARC00006');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00256E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00256E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000006',ix_);
        pb.xnames{iv} = 'PI000006';
        ig = ig_('ARC00007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.48655E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.48655E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000007',ix_);
        pb.xnames{iv} = 'PI000007';
        ig = ig_('ARC00007');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.48655E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.48655E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000007',ix_);
        pb.xnames{iv} = 'PI000007';
        ig = ig_('ARC00008');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.26895E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.26895E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000008',ix_);
        pb.xnames{iv} = 'PI000008';
        ig = ig_('ARC00010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.25622E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.25622E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000008',ix_);
        pb.xnames{iv} = 'PI000008';
        ig = ig_('ARC00011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.08033E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.08033E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000009',ix_);
        pb.xnames{iv} = 'PI000009';
        ig = ig_('ARC00010');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.25622E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.25622E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000009',ix_);
        pb.xnames{iv} = 'PI000009';
        ig = ig_('ARC00011');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.08033E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.08033E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000009',ix_);
        pb.xnames{iv} = 'PI000009';
        ig = ig_('ARC00012');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.81405E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.81405E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000009',ix_);
        pb.xnames{iv} = 'PI000009';
        ig = ig_('ARC00013');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.70084E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.70084E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000010',ix_);
        pb.xnames{iv} = 'PI000010';
        ig = ig_('ARC00012');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.81405E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.81405E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000010',ix_);
        pb.xnames{iv} = 'PI000010';
        ig = ig_('ARC00013');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.70084E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.70084E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000010',ix_);
        pb.xnames{iv} = 'PI000010';
        ig = ig_('ARC00014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.45124E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.45124E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000010',ix_);
        pb.xnames{iv} = 'PI000010';
        ig = ig_('ARC00015');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.16067E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.16067E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000011',ix_);
        pb.xnames{iv} = 'PI000011';
        ig = ig_('ARC00014');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.45124E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.45124E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000011',ix_);
        pb.xnames{iv} = 'PI000011';
        ig = ig_('ARC00015');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.16067E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.16067E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000011',ix_);
        pb.xnames{iv} = 'PI000011';
        ig = ig_('ARC00016');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.63836E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.63836E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000011',ix_);
        pb.xnames{iv} = 'PI000011';
        ig = ig_('ARC00021');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.14445E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.14445E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000012',ix_);
        pb.xnames{iv} = 'PI000012';
        ig = ig_('ARC00016');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.63836E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.63836E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000012',ix_);
        pb.xnames{iv} = 'PI000012';
        ig = ig_('ARC00017');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.07027E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.07027E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000013',ix_);
        pb.xnames{iv} = 'PI000013';
        ig = ig_('ARC00017');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.07027E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.07027E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000013',ix_);
        pb.xnames{iv} = 'PI000013';
        ig = ig_('ARC00018');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.25622E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.25622E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000014',ix_);
        pb.xnames{iv} = 'PI000014';
        ig = ig_('ARC00009');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.59656E-01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.59656E-01;
        end
        [iv,ix_] = s2mpjlib('ii','PI000014',ix_);
        pb.xnames{iv} = 'PI000014';
        ig = ig_('ARC00018');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 7.25622E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 7.25622E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000014',ix_);
        pb.xnames{iv} = 'PI000014';
        ig = ig_('ARC00019');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.62811E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.62811E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000015',ix_);
        pb.xnames{iv} = 'PI000015';
        ig = ig_('ARC00019');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.62811E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.62811E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000015',ix_);
        pb.xnames{iv} = 'PI000015';
        ig = ig_('ARC00020');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.45124E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.45124E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000016',ix_);
        pb.xnames{iv} = 'PI000016';
        ig = ig_('ARC00020');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.45124E+00+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.45124E+00;
        end
        [iv,ix_] = s2mpjlib('ii','PI000017',ix_);
        pb.xnames{iv} = 'PI000017';
        ig = ig_('ARC00021');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.14445E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.14445E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000017',ix_);
        pb.xnames{iv} = 'PI000017';
        ig = ig_('ARC00022');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.41977E-03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.41977E-03;
        end
        [iv,ix_] = s2mpjlib('ii','PI000018',ix_);
        pb.xnames{iv} = 'PI000018';
        ig = ig_('ARC00022');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.41977E-03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.41977E-03;
        end
        [iv,ix_] = s2mpjlib('ii','PI000018',ix_);
        pb.xnames{iv} = 'PI000018';
        ig = ig_('ARC00023');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.70320E-03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.70320E-03;
        end
        [iv,ix_] = s2mpjlib('ii','PI000019',ix_);
        pb.xnames{iv} = 'PI000019';
        ig = ig_('ARC00023');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.70320E-03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.70320E-03;
        end
        [iv,ix_] = s2mpjlib('ii','PI000019',ix_);
        pb.xnames{iv} = 'PI000019';
        ig = ig_('ARC00024');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.78190E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.78190E-02;
        end
        [iv,ix_] = s2mpjlib('ii','PI000020',ix_);
        pb.xnames{iv} = 'PI000020';
        ig = ig_('ARC00024');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.78190E-02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.78190E-02;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
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
        pbm.gconst(ig_('REGIO001')) = -3.91800E+00;
        pbm.gconst(ig_('REGIO002')) = -4.03400E+00;
        pbm.gconst(ig_('REGIO003')) = -5.25600E+00;
        pbm.gconst(ig_('REGIO004')) = -6.36500E+00;
        pbm.gconst(ig_('REGIO005')) = -2.12000E+00;
        pbm.gconst(ig_('REGIO006')) = -6.84800E+00;
        pbm.gconst(ig_('REGIO007')) = -1.56160E+01;
        pbm.gconst(ig_('REGIO008')) = -2.22000E-01;
        pbm.gconst(ig_('REGIO009')) = -1.91900E+00;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('FLOW0001'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0001')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0002'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0002')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0003'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0003')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0004'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0004')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0005'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0005')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0006'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0006')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0007'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0007')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0008'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0008')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0009'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0009')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0010'),1) = 0.00000E+00;
        pb.xupper(ix_('FLOW0010')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0011'),1) = 0.00000E+00;
        pb.xupper(ix_('FLOW0011')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0012'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0012')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0013'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0013')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0014'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0014')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0015'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0015')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0016'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0016')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0017'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0017')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0018'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0018')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0019'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0019')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0020'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0020')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0021'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0021')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0022'),1) = 0.00000E+00;
        pb.xupper(ix_('FLOW0022')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0023'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0023')) = 2.20120E+02;
        pb.xlower(ix_('FLOW0024'),1) = -2.20120E+02;
        pb.xupper(ix_('FLOW0024')) = 2.20120E+02;
        pb.xupper(ix_('SUPP0001')) = 1.15940E+01;
        pb.xupper(ix_('SUPP0002')) = 8.40000E+00;
        pb.xupper(ix_('SUPP0003')) = 4.80000E+00;
        pb.xupper(ix_('SUPP0004')) = 2.20120E+01;
        pb.xupper(ix_('SUPP0005')) = 1.20000E+00;
        pb.xupper(ix_('SUPP0006')) = 9.60000E-01;
        pb.xlower(ix_('PROD0001'),1) = 8.87000E+00;
        pb.xupper(ix_('PROD0001')) = 1.15940E+01;
        pb.xlower(ix_('PROD0002'),1) = 0.00000E+00;
        pb.xupper(ix_('PROD0002')) = 8.40000E+00;
        pb.xlower(ix_('PROD0003'),1) = 0.00000E+00;
        pb.xupper(ix_('PROD0003')) = 4.80000E+00;
        pb.xlower(ix_('PROD0004'),1) = 2.03440E+01;
        pb.xupper(ix_('PROD0004')) = 2.20120E+01;
        pb.xlower(ix_('PROD0005'),1) = 0.00000E+00;
        pb.xupper(ix_('PROD0005')) = 1.20000E+00;
        pb.xlower(ix_('PROD0006'),1) = 0.00000E+00;
        pb.xupper(ix_('PROD0006')) = 9.60000E-01;
        pb.xlower(ix_('PI000001'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000001')) = 5.92900E+03;
        pb.xlower(ix_('PI000002'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000002')) = 5.92900E+03;
        pb.xlower(ix_('PI000003'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000003')) = 6.40000E+03;
        pb.xlower(ix_('PI000004'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000004')) = 6.40000E+03;
        pb.xlower(ix_('PI000005'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000005')) = 5.92900E+03;
        pb.xlower(ix_('PI000006'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000006')) = 6.40000E+03;
        pb.xlower(ix_('PI000007'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000007')) = 6.40000E+03;
        pb.xlower(ix_('PI000008'),1) = 2.50000E+03;
        pb.xupper(ix_('PI000008')) = 4.38244E+03;
        pb.xlower(ix_('PI000009'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000009')) = 4.38244E+03;
        pb.xlower(ix_('PI000010'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000010')) = 4.38244E+03;
        pb.xlower(ix_('PI000011'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000011')) = 4.38244E+03;
        pb.xlower(ix_('PI000012'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000012')) = 4.38244E+03;
        pb.xlower(ix_('PI000013'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000013')) = 4.38244E+03;
        pb.xlower(ix_('PI000014'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000014')) = 4.38244E+03;
        pb.xlower(ix_('PI000015'),1) = 9.00000E+02;
        pb.xupper(ix_('PI000015')) = 4.38244E+03;
        pb.xlower(ix_('PI000016'),1) = 2.50000E+03;
        pb.xupper(ix_('PI000016')) = 4.38244E+03;
        pb.xlower(ix_('PI000017'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000017')) = 4.38244E+03;
        pb.xlower(ix_('PI000018'),1) = 0.00000E+00;
        pb.xupper(ix_('PI000018')) = 3.96900E+03;
        pb.xlower(ix_('PI000019'),1) = 6.25000E+02;
        pb.xupper(ix_('PI000019')) = 4.38244E+03;
        pb.xlower(ix_('PI000020'),1) = 6.25000E+02;
        pb.xupper(ix_('PI000020')) = 4.38244E+03;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('FLOW0001'),1) = 5.79700E+00;
        pb.x0(ix_('FLOW0002'),1) = 5.79700E+00;
        pb.x0(ix_('FLOW0003'),1) = 9.99700E+00;
        pb.x0(ix_('FLOW0004'),1) = 9.99700E+00;
        pb.x0(ix_('FLOW0005'),1) = 1.60760E+01;
        pb.x0(ix_('FLOW0006'),1) = 4.80000E+00;
        pb.x0(ix_('FLOW0007'),1) = 7.66000E-01;
        pb.x0(ix_('FLOW0008'),1) = -4.49000E+00;
        pb.x0(ix_('FLOW0009'),1) = 1.15860E+01;
        pb.x0(ix_('FLOW0010'),1) = 1.72404E+01;
        pb.x0(ix_('FLOW0011'),1) = 2.10363E+00;
        pb.x0(ix_('FLOW0012'),1) = 1.72404E+01;
        pb.x0(ix_('FLOW0013'),1) = 2.10363E+00;
        pb.x0(ix_('FLOW0014'),1) = 1.15676E+01;
        pb.x0(ix_('FLOW0015'),1) = 1.41145E+00;
        pb.x0(ix_('FLOW0016'),1) = 1.08380E+01;
        pb.x0(ix_('FLOW0017'),1) = 8.71800E+00;
        pb.x0(ix_('FLOW0018'),1) = 9.91800E+00;
        pb.x0(ix_('FLOW0019'),1) = 2.24640E+01;
        pb.x0(ix_('FLOW0020'),1) = 1.56160E+01;
        pb.x0(ix_('FLOW0021'),1) = 2.14100E+00;
        pb.x0(ix_('FLOW0022'),1) = 2.14100E+00;
        pb.x0(ix_('FLOW0023'),1) = 2.14100E+00;
        pb.x0(ix_('FLOW0024'),1) = 1.91900E+00;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'F00001SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0001';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00002SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0002';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00003SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0003';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00004SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0004';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00005SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0005';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00006SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0006';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00007SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0007';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00008SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0008';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00009SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0009';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00010SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0010';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00011SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0011';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00012SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0012';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00013SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0013';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00014SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0014';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00015SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0015';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00016SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0016';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00017SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0017';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00018SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0018';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00019SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0019';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00020SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0020';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00021SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0021';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00022SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0022';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00023SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0023';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'F00024SQ';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSQR';
        ielftype(ie) = iet_('eSQR');
        vname = 'FLOW0024';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('ARC00001');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00001SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00002');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00002SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00003');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00003SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00004');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00004SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00005');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00005SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00006');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00006SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00007');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00007SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00008');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00008SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00009');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00009SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00010');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00010SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00011');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00011SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00012');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00012SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00013');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00013SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00014');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00014SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00015');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00015SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00016');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00016SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00017');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00017SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00018');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00018SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00019');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00019SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00020');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00020SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00021');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00021SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00022');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00022SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00023');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00023SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('ARC00024');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('F00024SQ');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQI2-RN-65-59';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XGE0 = EV_(1)>=0.0e+0;
        if(XGE0)
            G = 2.0e+0*EV_(1);
        end
        if(~XGE0)
            G = -2.0e+0*EV_(1);
        end
        if(XGE0)
            H = 2.0e+0;
        end
        if(~XGE0)
            H = -2.0e+0;
        end
        varargout{1} = EV_(1)*abs(EV_(1));
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
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

