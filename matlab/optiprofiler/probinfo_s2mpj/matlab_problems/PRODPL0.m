function varargout = PRODPL0(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PRODPL0
%    *********
% 
%    A production planning problem in the computer industry.
% 
%    Source:
%    L. Escudero, private communication, 1991.
% 
%    SIF input: A.R. Conn, March 1991.
% 
%    classification = 'LQR2-RY-60-29'
% 
%    Constants
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PRODPL0';

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
        v_('T') = 5;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','COST',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','K01',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'K01';
        [ig,ig_] = s2mpjlib('ii','K02',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'K02';
        [ig,ig_] = s2mpjlib('ii','K03',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'K03';
        [ig,ig_] = s2mpjlib('ii','K04',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'K04';
        [ig,ig_] = s2mpjlib('ii','K05',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'K05';
        [ig,ig_] = s2mpjlib('ii','D00101',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00101';
        [ig,ig_] = s2mpjlib('ii','D00201',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00201';
        [ig,ig_] = s2mpjlib('ii','D00301',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00301';
        [ig,ig_] = s2mpjlib('ii','D00401',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00401';
        [ig,ig_] = s2mpjlib('ii','D00102',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00102';
        [ig,ig_] = s2mpjlib('ii','D00202',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00202';
        [ig,ig_] = s2mpjlib('ii','D00302',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00302';
        [ig,ig_] = s2mpjlib('ii','D00402',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00402';
        [ig,ig_] = s2mpjlib('ii','D00103',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00103';
        [ig,ig_] = s2mpjlib('ii','D00203',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00203';
        [ig,ig_] = s2mpjlib('ii','D00303',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00303';
        [ig,ig_] = s2mpjlib('ii','D00403',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00403';
        [ig,ig_] = s2mpjlib('ii','D00104',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00104';
        [ig,ig_] = s2mpjlib('ii','D00204',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00204';
        [ig,ig_] = s2mpjlib('ii','D00304',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00304';
        [ig,ig_] = s2mpjlib('ii','D00404',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00404';
        [ig,ig_] = s2mpjlib('ii','D00105',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00105';
        [ig,ig_] = s2mpjlib('ii','D00205',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00205';
        [ig,ig_] = s2mpjlib('ii','D00305',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00305';
        [ig,ig_] = s2mpjlib('ii','D00405',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'D00405';
        v_('TM1') = -1+v_('T');
        for I=v_('1'):v_('TM1')
            [ig,ig_] = s2mpjlib('ii',['SMOOTH',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['SMOOTH',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = numEntries(ig_);
        [iv,ix_] = s2mpjlib('ii','X00101',ix_);
        pb.xnames{iv} = 'X00101';
        ig = ig_('K01');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00201',ix_);
        pb.xnames{iv} = 'X00201';
        ig = ig_('K01');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00301',ix_);
        pb.xnames{iv} = 'X00301';
        ig = ig_('K01');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00401',ix_);
        pb.xnames{iv} = 'X00401';
        ig = ig_('K01');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00102',ix_);
        pb.xnames{iv} = 'X00102';
        ig = ig_('K02');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00202',ix_);
        pb.xnames{iv} = 'X00202';
        ig = ig_('K02');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00302',ix_);
        pb.xnames{iv} = 'X00302';
        ig = ig_('K02');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00402',ix_);
        pb.xnames{iv} = 'X00402';
        ig = ig_('K02');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00103',ix_);
        pb.xnames{iv} = 'X00103';
        ig = ig_('K03');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00203',ix_);
        pb.xnames{iv} = 'X00203';
        ig = ig_('K03');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00303',ix_);
        pb.xnames{iv} = 'X00303';
        ig = ig_('K03');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00403',ix_);
        pb.xnames{iv} = 'X00403';
        ig = ig_('K03');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00104',ix_);
        pb.xnames{iv} = 'X00104';
        ig = ig_('K04');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00204',ix_);
        pb.xnames{iv} = 'X00204';
        ig = ig_('K04');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00304',ix_);
        pb.xnames{iv} = 'X00304';
        ig = ig_('K04');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00404',ix_);
        pb.xnames{iv} = 'X00404';
        ig = ig_('K04');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00105',ix_);
        pb.xnames{iv} = 'X00105';
        ig = ig_('K05');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00205',ix_);
        pb.xnames{iv} = 'X00205';
        ig = ig_('K05');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00305',ix_);
        pb.xnames{iv} = 'X00305';
        ig = ig_('K05');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','X00405',ix_);
        pb.xnames{iv} = 'X00405';
        ig = ig_('K05');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        ig = ig_('D00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00101',ix_);
        pb.xnames{iv} = 'I00101';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000000;
        end
        ig = ig_('D00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00101',ix_);
        pb.xnames{iv} = 'I00101';
        ig = ig_('D00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00101',ix_);
        pb.xnames{iv} = 'Y00101';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00101');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00201',ix_);
        pb.xnames{iv} = 'I00201';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.000000;
        end
        ig = ig_('D00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00201',ix_);
        pb.xnames{iv} = 'I00201';
        ig = ig_('D00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00201',ix_);
        pb.xnames{iv} = 'Y00201';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3;
        end
        ig = ig_('D00201');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00301',ix_);
        pb.xnames{iv} = 'I00301';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.000000;
        end
        ig = ig_('D00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00301',ix_);
        pb.xnames{iv} = 'I00301';
        ig = ig_('D00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00301',ix_);
        pb.xnames{iv} = 'Y00301';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00301');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00401',ix_);
        pb.xnames{iv} = 'I00401';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.000000;
        end
        ig = ig_('D00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00401',ix_);
        pb.xnames{iv} = 'I00401';
        ig = ig_('D00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00401',ix_);
        pb.xnames{iv} = 'Y00401';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5;
        end
        ig = ig_('D00401');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00102',ix_);
        pb.xnames{iv} = 'I00102';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000000;
        end
        ig = ig_('D00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00102',ix_);
        pb.xnames{iv} = 'I00102';
        ig = ig_('D00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00102',ix_);
        pb.xnames{iv} = 'Y00102';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00102');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00202',ix_);
        pb.xnames{iv} = 'I00202';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.000000;
        end
        ig = ig_('D00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00202',ix_);
        pb.xnames{iv} = 'I00202';
        ig = ig_('D00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00202',ix_);
        pb.xnames{iv} = 'Y00202';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3;
        end
        ig = ig_('D00202');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00302',ix_);
        pb.xnames{iv} = 'I00302';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.000000;
        end
        ig = ig_('D00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00302',ix_);
        pb.xnames{iv} = 'I00302';
        ig = ig_('D00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00302',ix_);
        pb.xnames{iv} = 'Y00302';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00302');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00402',ix_);
        pb.xnames{iv} = 'I00402';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.000000;
        end
        ig = ig_('D00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00402',ix_);
        pb.xnames{iv} = 'I00402';
        ig = ig_('D00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00402',ix_);
        pb.xnames{iv} = 'Y00402';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5;
        end
        ig = ig_('D00402');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00103',ix_);
        pb.xnames{iv} = 'I00103';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000000;
        end
        ig = ig_('D00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00103',ix_);
        pb.xnames{iv} = 'I00103';
        ig = ig_('D00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00103',ix_);
        pb.xnames{iv} = 'Y00103';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00103');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00203',ix_);
        pb.xnames{iv} = 'I00203';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.000000;
        end
        ig = ig_('D00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00203',ix_);
        pb.xnames{iv} = 'I00203';
        ig = ig_('D00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00203',ix_);
        pb.xnames{iv} = 'Y00203';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3;
        end
        ig = ig_('D00203');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00303',ix_);
        pb.xnames{iv} = 'I00303';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.000000;
        end
        ig = ig_('D00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00303',ix_);
        pb.xnames{iv} = 'I00303';
        ig = ig_('D00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00303',ix_);
        pb.xnames{iv} = 'Y00303';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00303');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00403',ix_);
        pb.xnames{iv} = 'I00403';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.000000;
        end
        ig = ig_('D00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00403',ix_);
        pb.xnames{iv} = 'I00403';
        ig = ig_('D00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00403',ix_);
        pb.xnames{iv} = 'Y00403';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5;
        end
        ig = ig_('D00403');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00104',ix_);
        pb.xnames{iv} = 'I00104';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000000;
        end
        ig = ig_('D00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00104',ix_);
        pb.xnames{iv} = 'I00104';
        ig = ig_('D00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00104',ix_);
        pb.xnames{iv} = 'Y00104';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00104');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00204',ix_);
        pb.xnames{iv} = 'I00204';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.000000;
        end
        ig = ig_('D00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00204',ix_);
        pb.xnames{iv} = 'I00204';
        ig = ig_('D00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00204',ix_);
        pb.xnames{iv} = 'Y00204';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3;
        end
        ig = ig_('D00204');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00304',ix_);
        pb.xnames{iv} = 'I00304';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.000000;
        end
        ig = ig_('D00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00304',ix_);
        pb.xnames{iv} = 'I00304';
        ig = ig_('D00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00304',ix_);
        pb.xnames{iv} = 'Y00304';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00304');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00404',ix_);
        pb.xnames{iv} = 'I00404';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.000000;
        end
        ig = ig_('D00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00404',ix_);
        pb.xnames{iv} = 'I00404';
        ig = ig_('D00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00404',ix_);
        pb.xnames{iv} = 'Y00404';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5;
        end
        ig = ig_('D00404');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00105',ix_);
        pb.xnames{iv} = 'I00105';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.000000;
        end
        ig = ig_('D00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00105',ix_);
        pb.xnames{iv} = 'Y00105';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00105');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00205',ix_);
        pb.xnames{iv} = 'I00205';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.000000;
        end
        ig = ig_('D00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00205',ix_);
        pb.xnames{iv} = 'Y00205';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3;
        end
        ig = ig_('D00205');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00305',ix_);
        pb.xnames{iv} = 'I00305';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.000000;
        end
        ig = ig_('D00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00305',ix_);
        pb.xnames{iv} = 'Y00305';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2;
        end
        ig = ig_('D00305');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','I00405',ix_);
        pb.xnames{iv} = 'I00405';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.000000;
        end
        ig = ig_('D00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','Y00405',ix_);
        pb.xnames{iv} = 'Y00405';
        ig = ig_('COST');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5;
        end
        ig = ig_('D00405');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        pbm.gconst(ig_('K01')) = 3;
        pbm.gconst(ig_('K02')) = 6;
        pbm.gconst(ig_('K03')) = 10;
        pbm.gconst(ig_('K04')) = 2000;
        pbm.gconst(ig_('K05')) = 18;
        pbm.gconst(ig_('D00101')) = 1.000;
        pbm.gconst(ig_('D00201')) = 1.000;
        pbm.gconst(ig_('D00301')) = 1.000;
        pbm.gconst(ig_('D00401')) = 1.000;
        pbm.gconst(ig_('D00102')) = 2.667;
        pbm.gconst(ig_('D00202')) = 1.667;
        pbm.gconst(ig_('D00302')) = 2.667;
        pbm.gconst(ig_('D00402')) = 3.333;
        pbm.gconst(ig_('D00103')) = 2.667;
        pbm.gconst(ig_('D00203')) = 2.000;
        pbm.gconst(ig_('D00303')) = 3.000;
        pbm.gconst(ig_('D00403')) = 3.000;
        pbm.gconst(ig_('D00104')) = 2.667;
        pbm.gconst(ig_('D00204')) = 2.667;
        pbm.gconst(ig_('D00304')) = 2.667;
        pbm.gconst(ig_('D00404')) = 2.667;
        pbm.gconst(ig_('D00105')) = 2.667;
        pbm.gconst(ig_('D00205')) = 2.333;
        pbm.gconst(ig_('D00305')) = 2.333;
        pbm.gconst(ig_('D00405')) = 2.333;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQMRSQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('TM1')
            v_('IP1') = 1+I;
            ename = ['NLE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQMRSQ';
            ielftype(ie) = iet_('eSQMRSQ');
            vname = ['X0010',int2str(round(v_('IP1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0020',int2str(round(v_('IP1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0030',int2str(round(v_('IP1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0040',int2str(round(v_('IP1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0010',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0020',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0030',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X0040',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('TM1')
            ig = ig_(['SMOOTH',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['NLE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               58.7898356794
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-RY-60-29';
        pb.x0          = zeros(pb.n,1);
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 2.0;
        pbm.efpar(2) = 0.1;
        pbm.efpar(3) = pbm.efpar(2)*pbm.efpar(2);
        varargout{1} = pbm;

    case 'eSQMRSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,8);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        U_(1,3) = U_(1,3)+1;
        U_(1,4) = U_(1,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)+1;
        U_(2,7) = U_(2,7)+1;
        U_(2,8) = U_(2,8)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        U1MU2 = IV_(1)-IV_(2);
        varargout{1} = U1MU2^2-pbm.efpar(3)*IV_(2)^2;
        if(nargout>1)
            g_(1,1) = pbm.efpar(1)*U1MU2;
            g_(2,1) = -pbm.efpar(1)*U1MU2-pbm.efpar(1)*pbm.efpar(3)*IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.efpar(1);
                H_(1,2) = -pbm.efpar(1);
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.efpar(1)*(1.0-pbm.efpar(3));
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [3,0];
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

