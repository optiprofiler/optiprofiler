function varargout = NET1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NET1
%    *********
% 
%    A gas network problem for the south-east of England.
% 
%     SIF input: Sybille Schachler, Oxford, August 1992.
%    classification = 'OOI2-RN-48-57'
% 
%    ...Problem size parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NET1';

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
        v_('1') = 1;
        v_('NNOD') = 22;
        v_('NPIP') = 17;
        v_('NCMP') = 3;
        v_('NSRC') = 2;
        v_('CSTART') = 18;
        v_('CEND') = 20;
        v_('SSTART') = 21;
        v_('SEND') = 22;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','NOP17',ix_);
        pb.xnames{iv} = 'NOP17';
        [iv,ix_] = s2mpjlib('ii','PFL11',ix_);
        pb.xnames{iv} = 'PFL11';
        [iv,ix_] = s2mpjlib('ii','NOP9',ix_);
        pb.xnames{iv} = 'NOP9';
        [iv,ix_] = s2mpjlib('ii','PFL10',ix_);
        pb.xnames{iv} = 'PFL10';
        [iv,ix_] = s2mpjlib('ii','NOP16',ix_);
        pb.xnames{iv} = 'NOP16';
        [iv,ix_] = s2mpjlib('ii','PFL16',ix_);
        pb.xnames{iv} = 'PFL16';
        [iv,ix_] = s2mpjlib('ii','NOP19',ix_);
        pb.xnames{iv} = 'NOP19';
        [iv,ix_] = s2mpjlib('ii','PFL17',ix_);
        pb.xnames{iv} = 'PFL17';
        [iv,ix_] = s2mpjlib('ii','SFL22',ix_);
        pb.xnames{iv} = 'SFL22';
        [iv,ix_] = s2mpjlib('ii','SBV22',ix_);
        pb.xnames{iv} = 'SBV22';
        [iv,ix_] = s2mpjlib('ii','NOP18',ix_);
        pb.xnames{iv} = 'NOP18';
        [iv,ix_] = s2mpjlib('ii','NOP4',ix_);
        pb.xnames{iv} = 'NOP4';
        [iv,ix_] = s2mpjlib('ii','PFL13',ix_);
        pb.xnames{iv} = 'PFL13';
        [iv,ix_] = s2mpjlib('ii','PFL5',ix_);
        pb.xnames{iv} = 'PFL5';
        [iv,ix_] = s2mpjlib('ii','NOP11',ix_);
        pb.xnames{iv} = 'NOP11';
        [iv,ix_] = s2mpjlib('ii','PFL8',ix_);
        pb.xnames{iv} = 'PFL8';
        [iv,ix_] = s2mpjlib('ii','NOP6',ix_);
        pb.xnames{iv} = 'NOP6';
        [iv,ix_] = s2mpjlib('ii','NOP12',ix_);
        pb.xnames{iv} = 'NOP12';
        [iv,ix_] = s2mpjlib('ii','CFL19',ix_);
        pb.xnames{iv} = 'CFL19';
        [iv,ix_] = s2mpjlib('ii','CBV19',ix_);
        pb.xnames{iv} = 'CBV19';
        [iv,ix_] = s2mpjlib('ii','PFL7',ix_);
        pb.xnames{iv} = 'PFL7';
        [iv,ix_] = s2mpjlib('ii','NOP5',ix_);
        pb.xnames{iv} = 'NOP5';
        [iv,ix_] = s2mpjlib('ii','PFL6',ix_);
        pb.xnames{iv} = 'PFL6';
        [iv,ix_] = s2mpjlib('ii','NOP8',ix_);
        pb.xnames{iv} = 'NOP8';
        [iv,ix_] = s2mpjlib('ii','CFL20',ix_);
        pb.xnames{iv} = 'CFL20';
        [iv,ix_] = s2mpjlib('ii','CBV20',ix_);
        pb.xnames{iv} = 'CBV20';
        [iv,ix_] = s2mpjlib('ii','NOP7',ix_);
        pb.xnames{iv} = 'NOP7';
        [iv,ix_] = s2mpjlib('ii','PFL9',ix_);
        pb.xnames{iv} = 'PFL9';
        [iv,ix_] = s2mpjlib('ii','NOP21',ix_);
        pb.xnames{iv} = 'NOP21';
        [iv,ix_] = s2mpjlib('ii','PFL2',ix_);
        pb.xnames{iv} = 'PFL2';
        [iv,ix_] = s2mpjlib('ii','SFL21',ix_);
        pb.xnames{iv} = 'SFL21';
        [iv,ix_] = s2mpjlib('ii','SBV21',ix_);
        pb.xnames{iv} = 'SBV21';
        [iv,ix_] = s2mpjlib('ii','NOP1',ix_);
        pb.xnames{iv} = 'NOP1';
        [iv,ix_] = s2mpjlib('ii','PFL1',ix_);
        pb.xnames{iv} = 'PFL1';
        [iv,ix_] = s2mpjlib('ii','NOP14',ix_);
        pb.xnames{iv} = 'NOP14';
        [iv,ix_] = s2mpjlib('ii','PFL12',ix_);
        pb.xnames{iv} = 'PFL12';
        [iv,ix_] = s2mpjlib('ii','NOP10',ix_);
        pb.xnames{iv} = 'NOP10';
        [iv,ix_] = s2mpjlib('ii','PFL3',ix_);
        pb.xnames{iv} = 'PFL3';
        [iv,ix_] = s2mpjlib('ii','NOP2',ix_);
        pb.xnames{iv} = 'NOP2';
        [iv,ix_] = s2mpjlib('ii','CFL18',ix_);
        pb.xnames{iv} = 'CFL18';
        [iv,ix_] = s2mpjlib('ii','CBV18',ix_);
        pb.xnames{iv} = 'CBV18';
        [iv,ix_] = s2mpjlib('ii','NOP3',ix_);
        pb.xnames{iv} = 'NOP3';
        [iv,ix_] = s2mpjlib('ii','PFL4',ix_);
        pb.xnames{iv} = 'PFL4';
        [iv,ix_] = s2mpjlib('ii','NOP15',ix_);
        pb.xnames{iv} = 'NOP15';
        [iv,ix_] = s2mpjlib('ii','PFL15',ix_);
        pb.xnames{iv} = 'PFL15';
        [iv,ix_] = s2mpjlib('ii','NOP20',ix_);
        pb.xnames{iv} = 'NOP20';
        [iv,ix_] = s2mpjlib('ii','PFL14',ix_);
        pb.xnames{iv} = 'PFL14';
        [iv,ix_] = s2mpjlib('ii','NOP13',ix_);
        pb.xnames{iv} = 'NOP13';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','MBE1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE1';
        iv = ix_('PFL1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('PFL2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('SFL21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE2';
        iv = ix_('PFL3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CFL18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE3';
        iv = ix_('PFL4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CFL18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE4';
        iv = ix_('PFL5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('SFL22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE5';
        iv = ix_('PFL6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('PFL7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CFL19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE6';
        iv = ix_('PFL8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CFL19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE7';
        iv = ix_('PFL9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CFL20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE8';
        iv = ix_('PFL6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('CFL20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE9';
        iv = ix_('PFL10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('PFL11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE10';
        iv = ix_('PFL3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE11';
        iv = ix_('PFL5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE12';
        iv = ix_('PFL7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE13';
        iv = ix_('PFL14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE14';
        iv = ix_('PFL1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE15';
        iv = ix_('PFL4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE16';
        iv = ix_('PFL10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE17',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE17';
        iv = ix_('PFL11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE18';
        iv = ix_('PFL13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE19',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE19';
        iv = ix_('PFL16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE20';
        iv = ix_('PFL14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MBE21',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBE21';
        iv = ix_('PFL2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('PFL9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','MCR18',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MCR18';
        iv = ix_('NOP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000000;
        end
        iv = ix_('NOP2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.40000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.40000000;
        end
        iv = ix_('CBV18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00000;
        end
        [ig,ig_] = s2mpjlib('ii','MCR19',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MCR19';
        iv = ix_('NOP6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000000;
        end
        iv = ix_('NOP5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.40000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.40000000;
        end
        iv = ix_('CBV19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00000;
        end
        [ig,ig_] = s2mpjlib('ii','MCR20',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'MCR20';
        iv = ix_('NOP8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000000;
        end
        iv = ix_('NOP7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.40000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.40000000;
        end
        iv = ix_('CBV20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00000;
        end
        [ig,ig_] = s2mpjlib('ii','CLF18',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLF18';
        iv = ix_('CBV18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000e+04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000e+04;
        end
        iv = ix_('CFL18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','CLF19',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLF19';
        iv = ix_('CBV19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000e+04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000e+04;
        end
        iv = ix_('CFL19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','CLF20',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLF20';
        iv = ix_('CBV20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000e+04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000e+04;
        end
        iv = ix_('CFL20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','SLF21',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SLF21';
        iv = ix_('SBV21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00000;
        end
        iv = ix_('SFL21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','SUF21',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SUF21';
        iv = ix_('SBV21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.00000e+03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.00000e+03;
        end
        iv = ix_('SFL21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','SLF22',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SLF22';
        iv = ix_('SBV22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00000;
        end
        iv = ix_('SFL22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','SUF22',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'SUF22';
        iv = ix_('SBV22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.06000e+02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.06000e+02;
        end
        iv = ix_('SFL22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.00000000+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.00000000;
        end
        [ig,ig_] = s2mpjlib('ii','CLP18',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLP18';
        iv = ix_('NOP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('NOP2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CBV18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.72000e+02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.72000e+02;
        end
        [ig,ig_] = s2mpjlib('ii','CUP18',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CUP18';
        iv = ix_('NOP3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('NOP2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','CLP19',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLP19';
        iv = ix_('NOP6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('NOP5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CBV19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.45000e+02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.45000e+02;
        end
        [ig,ig_] = s2mpjlib('ii','CUP19',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CUP19';
        iv = ix_('NOP6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('NOP5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        [ig,ig_] = s2mpjlib('ii','CLP20',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CLP20';
        iv = ix_('NOP8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        iv = ix_('NOP7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('CBV20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.75000e+02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.75000e+02;
        end
        [ig,ig_] = s2mpjlib('ii','CUP20',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'CUP20';
        iv = ix_('NOP8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1;
        end
        iv = ix_('NOP7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1;
        end
        for i=v_('1'):v_('NPIP')
            [ig,ig_] = s2mpjlib('ii',['PDE',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PDE',int2str(i)];
            pbm.gscale(ig,1) = 1.00000e+03;
        end
        for i=v_('CSTART'):v_('CEND')
            [ig,ig_] = s2mpjlib('ii',['HPCON',int2str(i)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['HPCON',int2str(i)];
            pbm.gscale(ig,1) = 70.00000000;
        end
        for i=v_('CSTART'):v_('CEND')
            [ig,ig_] = s2mpjlib('ii',['HPOBJ',int2str(i)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = 0.03500000;
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
        pbm.gconst(ig_('MBE1')) = 6.65680000;
        pbm.gconst(ig_('MBE4')) = 1.96100000;
        pbm.gconst(ig_('MBE9')) = 3.72060e+02;
        pbm.gconst(ig_('MBE10')) = 47.17000000;
        pbm.gconst(ig_('MBE11')) = 1.60060e+02;
        pbm.gconst(ig_('MBE12')) = 4.25060e+02;
        pbm.gconst(ig_('MBE13')) = 5.30000e+02;
        pbm.gconst(ig_('MBE14')) = 24.16800000;
        pbm.gconst(ig_('MBE15')) = 2.54400000;
        pbm.gconst(ig_('MBE16')) = 89.14600000;
        pbm.gconst(ig_('MBE17')) = 4.92900e+02;
        pbm.gconst(ig_('MBE20')) = 4.64280e+02;
        pbm.gconst(ig_('MBE21')) = 1.48400e+02;
        pbm.gconst(ig_('CLF18')) = 1.00000e+04;
        pbm.gconst(ig_('CLF19')) = 1.00000e+04;
        pbm.gconst(ig_('CLF20')) = 1.00000e+04;
        pbm.gconst(ig_('HPCON18')) = 2.07000e+04;
        pbm.gconst(ig_('HPCON19')) = 2.07000e+04;
        pbm.gconst(ig_('HPCON20')) = 4.14000e+04;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('SFL21'),1) = 0.00000;
        pb.xupper(ix_('SFL21')) = 3.00000e+03;
        pb.xlower(ix_('SFL22'),1) = 0.00000;
        pb.xupper(ix_('SFL22')) = 1.06000e+02;
        pb.xlower(ix_('NOP1'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP1')) = 1.01500e+03;
        pb.xlower(ix_('NOP2'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP2')) = 1.10000e+03;
        pb.xlower(ix_('NOP3'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP3')) = 9.72000e+02;
        pb.xlower(ix_('NOP4'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP4')) = 1.10000e+03;
        pb.xlower(ix_('NOP5'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP5')) = 1.10000e+03;
        pb.xlower(ix_('NOP6'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP6')) = 8.45000e+02;
        pb.xlower(ix_('NOP7'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP7')) = 1.10000e+03;
        pb.xlower(ix_('NOP8'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP8')) = 1.07500e+03;
        pb.xlower(ix_('NOP9'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP9')) = 1.10000e+03;
        pb.xlower(ix_('NOP10'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP10')) = 1.10000e+03;
        pb.xlower(ix_('NOP11'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP11')) = 1.10000e+03;
        pb.xlower(ix_('NOP12'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP12')) = 1.10000e+03;
        pb.xlower(ix_('NOP13'),1) = 5.80000e+02;
        pb.xupper(ix_('NOP13')) = 1.10000e+03;
        pb.xlower(ix_('NOP14'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP14')) = 1.10000e+03;
        pb.xlower(ix_('NOP15'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP15')) = 1.10000e+03;
        pb.xlower(ix_('NOP16'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP16')) = 1.10000e+03;
        pb.xlower(ix_('NOP17'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP17')) = 1.10000e+03;
        pb.xlower(ix_('NOP18'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP18')) = 1.10000e+03;
        pb.xlower(ix_('NOP19'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP19')) = 1.10000e+03;
        pb.xlower(ix_('NOP20'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP20')) = 1.10000e+03;
        pb.xlower(ix_('NOP21'),1) = 5.00000e+02;
        pb.xupper(ix_('NOP21')) = 1.10000e+03;
        pb.xlower(ix_('CBV18'),1) = 0;
        pb.xupper(ix_('CBV18'),1) = 0;
        pb.xlower(ix_('CBV19'),1) = 1;
        pb.xupper(ix_('CBV19'),1) = 1;
        pb.xlower(ix_('CBV20'),1) = 1;
        pb.xupper(ix_('CBV20'),1) = 1;
        pb.xlower(ix_('SBV21'),1) = 1;
        pb.xupper(ix_('SBV21'),1) = 1;
        pb.xlower(ix_('SBV22'),1) = 1;
        pb.xupper(ix_('SBV22'),1) = 1;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 5.00000e+02*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eA0PANHAN',iet_);
        elftv{it}{1} = 'PIN';
        elftv{it}{2} = 'POUT';
        elftv{it}{3} = 'FLOW';
        elftp{it}{1} = 'PIPRES';
        [it,iet_] = s2mpjlib( 'ii', 'eA1MAXHP',iet_);
        elftv{it}{1} = 'PIN';
        elftv{it}{2} = 'POUT';
        elftv{it}{3} = 'FLOW';
        elftv{it}{4} = 'CBV';
        elftp{it}{1} = 'IPL';
        elftp{it}{2} = 'OPL';
        [it,iet_] = s2mpjlib( 'ii', 'eA2HPFUN',iet_);
        elftv{it}{1} = 'PIN';
        elftv{it}{2} = 'POUT';
        elftv{it}{3} = 'FLOW';
        elftv{it}{4} = 'CBV';
        elftp{it}{1} = 'IPL';
        elftp{it}{2} = 'OPL';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'PANH1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.62131268;
        ename = 'PANH2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.31605264;
        ename = 'PANH3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.13104611;
        ename = 'PANH4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.12796251;
        ename = 'PANH5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.78624623;
        ename = 'PANH6';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.84948702;
        ename = 'PANH7';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.13696026;
        ename = 'PANH8';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.25900862;
        ename = 'PANH9';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.32838618;
        ename = 'PANH10';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.33657520;
        ename = 'PANH11';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.61512113;
        ename = 'PANH12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.51339271;
        ename = 'PANH13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.20890923;
        ename = 'PANH14';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.15474706;
        ename = 'PANH15';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.26980036;
        ename = 'PANH16';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.04255562;
        ename = 'PANH17';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA0PANHAN';
        ielftype(ie) = iet_('eA0PANHAN');
        vname = 'NOP18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PFL17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PIPRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.12570329;
        ename = 'HPMAX18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA1MAXHP';
        ielftype(ie) = iet_('eA1MAXHP');
        vname = 'NOP2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        ename = 'HPMAX19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA1MAXHP';
        ielftype(ie) = iet_('eA1MAXHP');
        vname = 'NOP5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        ename = 'HPMAX20';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA1MAXHP';
        ielftype(ie) = iet_('eA1MAXHP');
        vname = 'NOP7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        ename = 'HPFUN18';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA2HPFUN';
        ielftype(ie) = iet_('eA2HPFUN');
        vname = 'NOP2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        ename = 'HPFUN19';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA2HPFUN';
        ielftype(ie) = iet_('eA2HPFUN');
        vname = 'NOP5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        ename = 'HPFUN20';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eA2HPFUN';
        ielftype(ie) = iet_('eA2HPFUN');
        vname = 'NOP7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('PIN',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'NOP8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('POUT',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CFL20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('FLOW',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'CBV20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],5.00000e+02);
        posev = find(strcmp('CBV',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('IPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        [~,posep] = ismember('OPL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00000;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for i=v_('1'):v_('NPIP')
            ig = ig_(['PDE',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PANH',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for i=v_('CSTART'):v_('CEND')
            ig = ig_(['HPCON',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['HPMAX',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for i=v_('CSTART'):v_('CEND')
            ig = ig_(['HPOBJ',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['HPFUN',int2str(i)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOI2-RN-48-57';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eA0PANHAN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A0FLEX = 1.8539e0;
        A0FGT0 = EV_(3)>=0.0e0;
        if(A0FGT0)
            A0HFLO =...
              -pbm.elpar{iel_}(1)*A0FLEX*(A0FLEX-1.0e0)*EV_(3)^(A0FLEX-2.0e0);
        end
        if(~A0FGT0)
            A0HFLO =...
              A0FLEX*(A0FLEX-1.0e0)*pbm.elpar{iel_}(1)*abs(EV_(3))^(A0FLEX-2.0e0);
        end
        varargout{1} = EV_(1)*EV_(1)-EV_(2)*EV_(2)-pbm.elpar{iel_}(1)*EV_(3)*...
             abs(EV_(3))^(A0FLEX-1.0e0);
        if(nargout>1)
            g_(2,1) = -2.0e0*EV_(2);
            g_(1,1) = 2.0e0*EV_(1);
            g_(3,1) = -pbm.elpar{iel_}(1)*A0FLEX*abs(EV_(3))^(A0FLEX-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(2,2) = -2.0e0;
                H_(1,1) = 2.0e0;
                H_(3,3) = A0HFLO;
                varargout{3} = H_;
            end
        end

    case 'eA1MAXHP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A1BETA = 0.23077e0;
        A1HFAC = 203.712e0;
        A1PSUC = EV_(1)-pbm.elpar{iel_}(1)*EV_(4);
        A1PDIS = EV_(2)+pbm.elpar{iel_}(2)*EV_(4);
        A1CRB = (A1PDIS/A1PSUC)^A1BETA;
        A1PROD = A1BETA*A1HFAC*A1CRB*EV_(3);
        A1GPIN = -A1PROD/A1PSUC;
        A1GPOU = A1PROD/A1PDIS;
        A1GFLO = A1HFAC*(A1CRB-1.0e0);
        A1GCBV = -pbm.elpar{iel_}(1)*A1GPIN+pbm.elpar{iel_}(2)*A1GPOU;
        A1HII = A1PROD*(A1BETA+1.0e0)/(A1PSUC^2);
        A1HIO = -A1PROD*A1BETA/(A1PSUC*A1PDIS);
        A1HOO = A1PROD*(A1BETA-1.0e0)/(A1PDIS^2);
        A1HIC = -pbm.elpar{iel_}(1)*A1HII+pbm.elpar{iel_}(2)*A1HIO;
        A1HOC = -pbm.elpar{iel_}(1)*A1HIO+pbm.elpar{iel_}(2)*A1HOO;
        varargout{1} = EV_(3)*A1GFLO;
        if(nargout>1)
            g_(1,1) = A1GPIN;
            g_(2,1) = A1GPOU;
            g_(3,1) = A1GFLO;
            g_(4,1) = A1GCBV;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = A1HII;
                H_(1,2) = A1HIO;
                H_(2,1) = H_(1,2);
                H_(2,2) = A1HOO;
                H_(1,3) = A1GPIN/EV_(3);
                H_(3,1) = H_(1,3);
                H_(2,3) = A1GPOU/EV_(3);
                H_(3,2) = H_(2,3);
                H_(1,4) = A1HIC;
                H_(4,1) = H_(1,4);
                H_(2,4) = A1HOC;
                H_(4,2) = H_(2,4);
                H_(3,4) = A1GCBV/EV_(3);
                H_(4,3) = H_(3,4);
                H_(4,4) = -pbm.elpar{iel_}(1)*A1HIC+pbm.elpar{iel_}(2)*A1HOC;
                varargout{3} = H_;
            end
        end

    case 'eA2HPFUN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A2BETA = 0.23077e0;
        A2HFAC = 203.712e0;
        A2PSUC = EV_(1)-pbm.elpar{iel_}(1)*EV_(4);
        A2PDIS = EV_(2)+pbm.elpar{iel_}(2)*EV_(4);
        A2CRB = (A2PDIS/A2PSUC)^A2BETA;
        A2PROD = A2BETA*A2HFAC*A2CRB*EV_(3);
        A2GPIN = -A2PROD/A2PSUC;
        A2GPOU = A2PROD/A2PDIS;
        A2GFLO = A2HFAC*(A2CRB-1.0e0);
        A2GCBV = -pbm.elpar{iel_}(1)*A2GPIN+pbm.elpar{iel_}(2)*A2GPOU;
        A2HII = A2PROD*(A2BETA+1.0e0)/(A2PSUC^2);
        A2HIO = -A2PROD*A2BETA/(A2PSUC*A2PDIS);
        A2HOO = A2PROD*(A2BETA-1.0e0)/(A2PDIS^2);
        A2HIC = -pbm.elpar{iel_}(1)*A2HII+pbm.elpar{iel_}(2)*A2HIO;
        A2HOC = -pbm.elpar{iel_}(1)*A2HIO+pbm.elpar{iel_}(2)*A2HOO;
        varargout{1} = EV_(3)*A2GFLO;
        if(nargout>1)
            g_(1,1) = A2GPIN;
            g_(2,1) = A2GPOU;
            g_(3,1) = A2GFLO;
            g_(4,1) = A2GCBV;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = A2HII;
                H_(1,2) = A2HIO;
                H_(2,1) = H_(1,2);
                H_(2,2) = A2HOO;
                H_(1,3) = A2GPIN/EV_(3);
                H_(3,1) = H_(1,3);
                H_(2,3) = A2GPOU/EV_(3);
                H_(3,2) = H_(2,3);
                H_(1,4) = A2HIC;
                H_(4,1) = H_(1,4);
                H_(2,4) = A2HOC;
                H_(4,2) = H_(2,4);
                H_(3,4) = A2GCBV/EV_(3);
                H_(4,3) = H_(3,4);
                H_(4,4) = -pbm.elpar{iel_}(1)*A2HIC+pbm.elpar{iel_}(2)*A2HOC;
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

