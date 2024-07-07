function varargout = DEGENLPA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DEGENLPA
%    *********
% 
%    A small linear program with some degeneracy.
% 
%    Source:
%    T.C.T. Kotiah and D.I. Steinberg,
%    "Occurences of cycling and other phenomena arising in a class of
%    linear programming models",
%    Communications of the ACM, vol. 20, pp. 107-112, 1977.
% 
%    SIF input: Ph. Toint, Aug 1990.
% 
%    classification = 'LLR2-AN-20-15'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DEGENLPA';

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
        v_('N') = 20;
        v_('M') = 15;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 33.333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 33.333;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 33.343+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 33.343;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.01;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 33.333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 33.333;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 133.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 133.33;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        [ig,ig_] = s2mpjlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 300.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 300.0;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.09+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.09;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        [ig,ig_] = s2mpjlib('ii','C3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C3';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.326+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.326;
        end
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -101.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -101.0;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 200.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 200.0;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.06;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02;
        end
        [ig,ig_] = s2mpjlib('ii','C4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C4';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.0066667+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.0066667;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.03;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 200.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 200.0;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.06;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02;
        end
        [ig,ig_] = s2mpjlib('ii','C5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C5';
        iv = ix_('X1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.6667e-4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.6667e-4;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.01;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 200.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 200.0;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.06;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02;
        end
        [ig,ig_] = s2mpjlib('ii','C6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C6';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.978+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.978;
        end
        iv = ix_('X5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -201.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -201.0;
        end
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C7';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.489+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.489;
        end
        iv = ix_('X6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -101.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -101.03;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C8';
        iv = ix_('X2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.001;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.489+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.489;
        end
        iv = ix_('X7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -101.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -101.03;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C9';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.001+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.001;
        end
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        iv = ix_('X9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.04;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C10';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.02;
        end
        iv = ix_('X8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.06;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C11';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.002+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.002;
        end
        iv = ix_('X10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.02+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.02;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 100.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 100.0;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.03;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.01;
        end
        [ig,ig_] = s2mpjlib('ii','C12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C12';
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.5742e-6+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.5742e-6;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00252+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00252;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.61975+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.61975;
        end
        iv = ix_('X20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.03;
        end
        [ig,ig_] = s2mpjlib('ii','C13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C13';
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.00257+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.00257;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.25221+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.25221;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -6.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -6.2;
        end
        iv = ix_('X17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.09+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.09;
        end
        [ig,ig_] = s2mpjlib('ii','C14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C14';
        iv = ix_('X11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00629+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00629;
        end
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.20555+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.20555;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.1106+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.1106;
        end
        iv = ix_('X15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 101.04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 101.04;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 505.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 505.1;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -256.72+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -256.72;
        end
        [ig,ig_] = s2mpjlib('ii','C15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C15';
        iv = ix_('X12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.00841+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.00841;
        end
        iv = ix_('X13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.08406+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.08406;
        end
        iv = ix_('X14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -0.20667+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -0.20667;
        end
        iv = ix_('X16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 20.658+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 20.658;
        end
        iv = ix_('X18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.07+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.07;
        end
        iv = ix_('X19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.5;
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
        pbm.gconst(ig_('C1')) = 0.70785;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 1.0*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               3.06435
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-AN-20-15';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

