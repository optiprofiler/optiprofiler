function varargout = QPCBLEND(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : QPCBLEND
%    *********
% 
%    Source: a variant on the BLEND linear programming problem
%    with an additional CONVEX diagonal Hessian matrix as given by
%    N. I. M. Gould, "An algorithm for large-scale quadratic programming",
%    IMA J. Num. Anal (1991), 11, 299-324, problem class 4.
% 
%    SIF input: Nick Gould, January 1993
% 
%    classification = 'C-CQLR2-MN-83-74'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'QPCBLEND';

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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '1';
        [ig,ig_] = s2mpjlib('ii','2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '2';
        [ig,ig_] = s2mpjlib('ii','3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '3';
        [ig,ig_] = s2mpjlib('ii','4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '4';
        [ig,ig_] = s2mpjlib('ii','5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '5';
        [ig,ig_] = s2mpjlib('ii','6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '6';
        [ig,ig_] = s2mpjlib('ii','7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '7';
        [ig,ig_] = s2mpjlib('ii','8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '8';
        [ig,ig_] = s2mpjlib('ii','9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '9';
        [ig,ig_] = s2mpjlib('ii','10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '10';
        [ig,ig_] = s2mpjlib('ii','11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '11';
        [ig,ig_] = s2mpjlib('ii','12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '12';
        [ig,ig_] = s2mpjlib('ii','13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '13';
        [ig,ig_] = s2mpjlib('ii','14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '14';
        [ig,ig_] = s2mpjlib('ii','15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '15';
        [ig,ig_] = s2mpjlib('ii','16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '16';
        [ig,ig_] = s2mpjlib('ii','17',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '17';
        [ig,ig_] = s2mpjlib('ii','18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '18';
        [ig,ig_] = s2mpjlib('ii','19',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '19';
        [ig,ig_] = s2mpjlib('ii','20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '20';
        [ig,ig_] = s2mpjlib('ii','21',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '21';
        [ig,ig_] = s2mpjlib('ii','22',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '22';
        [ig,ig_] = s2mpjlib('ii','23',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '23';
        [ig,ig_] = s2mpjlib('ii','24',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '24';
        [ig,ig_] = s2mpjlib('ii','25',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '25';
        [ig,ig_] = s2mpjlib('ii','26',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '26';
        [ig,ig_] = s2mpjlib('ii','27',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '27';
        [ig,ig_] = s2mpjlib('ii','28',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '28';
        [ig,ig_] = s2mpjlib('ii','29',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '29';
        [ig,ig_] = s2mpjlib('ii','30',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '30';
        [ig,ig_] = s2mpjlib('ii','31',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '31';
        [ig,ig_] = s2mpjlib('ii','32',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '32';
        [ig,ig_] = s2mpjlib('ii','33',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '33';
        [ig,ig_] = s2mpjlib('ii','34',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '34';
        [ig,ig_] = s2mpjlib('ii','35',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '35';
        [ig,ig_] = s2mpjlib('ii','36',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '36';
        [ig,ig_] = s2mpjlib('ii','37',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '37';
        [ig,ig_] = s2mpjlib('ii','38',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '38';
        [ig,ig_] = s2mpjlib('ii','39',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '39';
        [ig,ig_] = s2mpjlib('ii','40',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '40';
        [ig,ig_] = s2mpjlib('ii','41',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '41';
        [ig,ig_] = s2mpjlib('ii','42',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '42';
        [ig,ig_] = s2mpjlib('ii','43',ig_);
        gtype{ig}  = '==';
        cnames{ig} = '43';
        [ig,ig_] = s2mpjlib('ii','44',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '44';
        [ig,ig_] = s2mpjlib('ii','45',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '45';
        [ig,ig_] = s2mpjlib('ii','46',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '46';
        [ig,ig_] = s2mpjlib('ii','47',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '47';
        [ig,ig_] = s2mpjlib('ii','48',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '48';
        [ig,ig_] = s2mpjlib('ii','49',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '49';
        [ig,ig_] = s2mpjlib('ii','50',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '50';
        [ig,ig_] = s2mpjlib('ii','51',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '51';
        [ig,ig_] = s2mpjlib('ii','52',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '52';
        [ig,ig_] = s2mpjlib('ii','53',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '53';
        [ig,ig_] = s2mpjlib('ii','54',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '54';
        [ig,ig_] = s2mpjlib('ii','55',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '55';
        [ig,ig_] = s2mpjlib('ii','56',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '56';
        [ig,ig_] = s2mpjlib('ii','57',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '57';
        [ig,ig_] = s2mpjlib('ii','58',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '58';
        [ig,ig_] = s2mpjlib('ii','59',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '59';
        [ig,ig_] = s2mpjlib('ii','60',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '60';
        [ig,ig_] = s2mpjlib('ii','61',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '61';
        [ig,ig_] = s2mpjlib('ii','62',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '62';
        [ig,ig_] = s2mpjlib('ii','63',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '63';
        [ig,ig_] = s2mpjlib('ii','64',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '64';
        [ig,ig_] = s2mpjlib('ii','65',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '65';
        [ig,ig_] = s2mpjlib('ii','66',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '66';
        [ig,ig_] = s2mpjlib('ii','67',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '67';
        [ig,ig_] = s2mpjlib('ii','68',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '68';
        [ig,ig_] = s2mpjlib('ii','69',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '69';
        [ig,ig_] = s2mpjlib('ii','70',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '70';
        [ig,ig_] = s2mpjlib('ii','71',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '71';
        [ig,ig_] = s2mpjlib('ii','72',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '72';
        [ig,ig_] = s2mpjlib('ii','73',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '73';
        [ig,ig_] = s2mpjlib('ii','74',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = '74';
        [ig,ig_] = s2mpjlib('ii','C',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = ig_.Count;
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.537+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.537;
        end
        ig = ig_('3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.131+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.131;
        end
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1155;
        end
        ig = ig_('5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0365+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0365;
        end
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.143+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.143;
        end
        ig = ig_('7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.037+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.037;
        end
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .003;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0587+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0587;
        end
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .15+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .15;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .302+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .302;
        end
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        ig = ig_('67');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.2;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2931+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2931;
        end
        ig = ig_('3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.117+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.117;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0649+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0649;
        end
        ig = ig_('5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1233;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2217+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2217;
        end
        ig = ig_('8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.18+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.18;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0042+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0042;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .003+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .003;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1053+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1053;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .185;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .384+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .384;
        end
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.00862+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.00862;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.00862+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.00862;
        end
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0101+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0101;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0101+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0101;
        end
        ig = ig_('68');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.87+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.87;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0277+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0277;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0563+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0563;
        end
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.199+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.199;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.6873+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.6873;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.017+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.017;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01303;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0506+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0506;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .209+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .209;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .495+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .495;
        end
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0112+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0112;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0378+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0378;
        end
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1502+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1502;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.7953+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.7953;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0099+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0099;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01303;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0448+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0448;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .185;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .721+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .721;
        end
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.175+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.175;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.27+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.27;
        end
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.028+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.028;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.455+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.455;
        end
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        ig = ig_('21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01303;
        end
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0506+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0506;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .209+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .209;
        end
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .495+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .495;
        end
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.271+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.271;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.3285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.3285;
        end
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0255+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0255;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2656+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2656;
        end
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        ig = ig_('18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01303;
        end
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0506+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0506;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .209+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .209;
        end
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .495+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .495;
        end
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2836+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2836;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.3285+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.3285;
        end
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0241;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2502+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2502;
        end
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        ig = ig_('17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01303;
        end
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0506+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0506;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .209+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .209;
        end
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .495+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .495;
        end
        [iv,ix_] = s2mpjlib('ii','8',ix_);
        pb.xnames{iv} = '8';
        ig = ig_('12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','8',ix_);
        pb.xnames{iv} = '8';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0327+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0327;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .094+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .094;
        end
        [iv,ix_] = s2mpjlib('ii','8',ix_);
        pb.xnames{iv} = '8';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .045+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .045;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .793+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .793;
        end
        [iv,ix_] = s2mpjlib('ii','8',ix_);
        pb.xnames{iv} = '8';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0044+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0044;
        end
        [iv,ix_] = s2mpjlib('ii','9',ix_);
        pb.xnames{iv} = '9';
        ig = ig_('15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','10',ix_);
        pb.xnames{iv} = '10';
        ig = ig_('16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','11',ix_);
        pb.xnames{iv} = '11';
        ig = ig_('14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','12',ix_);
        pb.xnames{iv} = '12';
        ig = ig_('14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('15');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('17');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0588+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0588;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.8145+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.8145;
        end
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0091+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0091;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.8239+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.8239;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0081+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0081;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2112+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2112;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .387+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .387;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.03;
        end
        ig = ig_('69');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.3;
        end
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .07+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .07;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('16');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('18');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0404+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0404;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.8564+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.8564;
        end
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0069+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0069;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.7689+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.7689;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0063+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0063;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.156+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.156;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .297+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .297;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .792+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .792;
        end
        ig = ig_('69');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0378+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0378;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.3321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.3321;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5875+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5875;
        end
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.362;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.3;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2049+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2049;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .826;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.61+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.61;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('70');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .155;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.3321+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.3321;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5875+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5875;
        end
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.362+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.362;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.3;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2049+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2049;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .826;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.61+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.61;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('66');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('70');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .155;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2414+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2414;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.6627+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.6627;
        end
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.293+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.293;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.3;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1531+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1531;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .826;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.61+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.61;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('70');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .155;
        end
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        ig = ig_('21');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2414+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2414;
        end
        ig = ig_('22');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.6627+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.6627;
        end
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.293+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.293;
        end
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.3;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.1531+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.1531;
        end
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .826+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .826;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.61+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.61;
        end
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        ig = ig_('70');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .155+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .155;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0185;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0568+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0568;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0806+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0806;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0658+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0658;
        end
        ig = ig_('26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0328+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0328;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.4934+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.4934;
        end
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2922+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2922;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0096+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0096;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0654;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2535;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .632;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .6807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .6807;
        end
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        ig = ig_('71');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0528+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0528;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0185+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0185;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0568+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0568;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0806+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0806;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0658+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0658;
        end
        ig = ig_('26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0328+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0328;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.4934+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.4934;
        end
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2922+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2922;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0096+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0096;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0654;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2535+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2535;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .632;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .6807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .6807;
        end
        ig = ig_('66');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        ig = ig_('71');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0528+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0528;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0184+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0184;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0564+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0564;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.078;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0655+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0655;
        end
        ig = ig_('26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0303;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.475;
        end
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.305+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.305;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0654;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2703+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2703;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .632;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .6807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .6807;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('71');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0528+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0528;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0184+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0184;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0564+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0564;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.078+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.078;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0655+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0655;
        end
        ig = ig_('26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0303+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0303;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.475+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.475;
        end
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.305+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.305;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0654+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0654;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.2703+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.2703;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .632+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .632;
        end
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .6807+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .6807;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('71');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0528+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0528;
        end
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .76+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .76;
        end
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .5714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .5714;
        end
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        ig = ig_('30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1869+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1869;
        end
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .2796+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .2796;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.241+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.241;
        end
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.766+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.766;
        end
        ig = ig_('72');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .128+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .128;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0571+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0571;
        end
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.0114+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.0114;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .6571+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .6571;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .5714+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .5714;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .1724+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .1724;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .2579+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .2579;
        end
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.067+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.067;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.552+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.552;
        end
        ig = ig_('72');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .118+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .118;
        end
        [iv,ix_] = s2mpjlib('ii','25',ix_);
        pb.xnames{iv} = '25';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','26',ix_);
        pb.xnames{iv} = '26';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','27',ix_);
        pb.xnames{iv} = '27';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','28',ix_);
        pb.xnames{iv} = '28';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','28',ix_);
        pb.xnames{iv} = '28';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.95+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.95;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.7+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.7;
        end
        [iv,ix_] = s2mpjlib('ii','28',ix_);
        pb.xnames{iv} = '28';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.0;
        end
        [iv,ix_] = s2mpjlib('ii','28',ix_);
        pb.xnames{iv} = '28';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','29',ix_);
        pb.xnames{iv} = '29';
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','29',ix_);
        pb.xnames{iv} = '29';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.84+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.84;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.45+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.45;
        end
        [iv,ix_] = s2mpjlib('ii','29',ix_);
        pb.xnames{iv} = '29';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 12.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 12.0;
        end
        [iv,ix_] = s2mpjlib('ii','29',ix_);
        pb.xnames{iv} = '29';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','30',ix_);
        pb.xnames{iv} = '30';
        ig = ig_('19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','30',ix_);
        pb.xnames{iv} = '30';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.43+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.43;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.57+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.57;
        end
        [iv,ix_] = s2mpjlib('ii','30',ix_);
        pb.xnames{iv} = '30';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','30',ix_);
        pb.xnames{iv} = '30';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .233;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.358+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.358;
        end
        [iv,ix_] = s2mpjlib('ii','31',ix_);
        pb.xnames{iv} = '31';
        ig = ig_('20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','31',ix_);
        pb.xnames{iv} = '31';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.03;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.32+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.32;
        end
        [iv,ix_] = s2mpjlib('ii','31',ix_);
        pb.xnames{iv} = '31';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','31',ix_);
        pb.xnames{iv} = '31';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .205;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.333;
        end
        [iv,ix_] = s2mpjlib('ii','32',ix_);
        pb.xnames{iv} = '32';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','32',ix_);
        pb.xnames{iv} = '32';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.23+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.23;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.22+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.22;
        end
        [iv,ix_] = s2mpjlib('ii','32',ix_);
        pb.xnames{iv} = '32';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        [iv,ix_] = s2mpjlib('ii','32',ix_);
        pb.xnames{iv} = '32';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .381+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .381;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.509+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.509;
        end
        [iv,ix_] = s2mpjlib('ii','33',ix_);
        pb.xnames{iv} = '33';
        ig = ig_('30');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','33',ix_);
        pb.xnames{iv} = '33';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.4;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.85+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.85;
        end
        [iv,ix_] = s2mpjlib('ii','33',ix_);
        pb.xnames{iv} = '33';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.5;
        end
        [iv,ix_] = s2mpjlib('ii','33',ix_);
        pb.xnames{iv} = '33';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .39+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .39;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.77+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.77;
        end
        [iv,ix_] = s2mpjlib('ii','34',ix_);
        pb.xnames{iv} = '34';
        ig = ig_('31');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','34',ix_);
        pb.xnames{iv} = '34';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.74+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.74;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.1;
        end
        [iv,ix_] = s2mpjlib('ii','34',ix_);
        pb.xnames{iv} = '34';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.3+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.3;
        end
        [iv,ix_] = s2mpjlib('ii','34',ix_);
        pb.xnames{iv} = '34';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .233;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.58+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.58;
        end
        [iv,ix_] = s2mpjlib('ii','35',ix_);
        pb.xnames{iv} = '35';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','35',ix_);
        pb.xnames{iv} = '35';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.74+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.74;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.9;
        end
        [iv,ix_] = s2mpjlib('ii','35',ix_);
        pb.xnames{iv} = '35';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.0;
        end
        [iv,ix_] = s2mpjlib('ii','35',ix_);
        pb.xnames{iv} = '35';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','36',ix_);
        pb.xnames{iv} = '36';
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.493+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.493;
        end
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.165+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.165;
        end
        [iv,ix_] = s2mpjlib('ii','36',ix_);
        pb.xnames{iv} = '36';
        ig = ig_('46');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0924+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0924;
        end
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        ig = ig_('32');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('44');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.03;
        end
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        ig = ig_('45');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.03+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.03;
        end
        ig = ig_('47');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.5;
        end
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        ig = ig_('48');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5;
        end
        ig = ig_('49');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .5;
        end
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        ig = ig_('73');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .64+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .64;
        end
        ig = ig_('74');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .35+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .35;
        end
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.36+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.36;
        end
        [iv,ix_] = s2mpjlib('ii','38',ix_);
        pb.xnames{iv} = '38';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','38',ix_);
        pb.xnames{iv} = '38';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.98+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.98;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.58+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.58;
        end
        [iv,ix_] = s2mpjlib('ii','38',ix_);
        pb.xnames{iv} = '38';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.0;
        end
        [iv,ix_] = s2mpjlib('ii','38',ix_);
        pb.xnames{iv} = '38';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','39',ix_);
        pb.xnames{iv} = '39';
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','39',ix_);
        pb.xnames{iv} = '39';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.87+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.87;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.33+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.33;
        end
        [iv,ix_] = s2mpjlib('ii','39',ix_);
        pb.xnames{iv} = '39';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 12.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 12.0;
        end
        [iv,ix_] = s2mpjlib('ii','39',ix_);
        pb.xnames{iv} = '39';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','40',ix_);
        pb.xnames{iv} = '40';
        ig = ig_('19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','40',ix_);
        pb.xnames{iv} = '40';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.46;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.45+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.45;
        end
        [iv,ix_] = s2mpjlib('ii','40',ix_);
        pb.xnames{iv} = '40';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','40',ix_);
        pb.xnames{iv} = '40';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .233;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.358+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.358;
        end
        [iv,ix_] = s2mpjlib('ii','41',ix_);
        pb.xnames{iv} = '41';
        ig = ig_('20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','41',ix_);
        pb.xnames{iv} = '41';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.06+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.06;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.2;
        end
        [iv,ix_] = s2mpjlib('ii','41',ix_);
        pb.xnames{iv} = '41';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','41',ix_);
        pb.xnames{iv} = '41';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .205;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.333;
        end
        [iv,ix_] = s2mpjlib('ii','42',ix_);
        pb.xnames{iv} = '42';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','42',ix_);
        pb.xnames{iv} = '42';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.26+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.26;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.13+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.13;
        end
        [iv,ix_] = s2mpjlib('ii','42',ix_);
        pb.xnames{iv} = '42';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        [iv,ix_] = s2mpjlib('ii','42',ix_);
        pb.xnames{iv} = '42';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .318+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .318;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.509+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.509;
        end
        [iv,ix_] = s2mpjlib('ii','43',ix_);
        pb.xnames{iv} = '43';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','43',ix_);
        pb.xnames{iv} = '43';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.77+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.77;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.78+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.78;
        end
        [iv,ix_] = s2mpjlib('ii','43',ix_);
        pb.xnames{iv} = '43';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.0;
        end
        [iv,ix_] = s2mpjlib('ii','43',ix_);
        pb.xnames{iv} = '43';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','44',ix_);
        pb.xnames{iv} = '44';
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.435+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.435;
        end
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.208+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.208;
        end
        [iv,ix_] = s2mpjlib('ii','44',ix_);
        pb.xnames{iv} = '44';
        ig = ig_('52');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0924+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0924;
        end
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        ig = ig_('33');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('50');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.65+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.65;
        end
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        ig = ig_('51');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.65+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.65;
        end
        ig = ig_('53');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.5;
        end
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        ig = ig_('54');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5;
        end
        ig = ig_('55');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .5;
        end
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        ig = ig_('73');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.36+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.36;
        end
        ig = ig_('74');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .35+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .35;
        end
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -5.08+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -5.08;
        end
        [iv,ix_] = s2mpjlib('ii','46',ix_);
        pb.xnames{iv} = '46';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','46',ix_);
        pb.xnames{iv} = '46';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -7.99+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -7.99;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.59+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.59;
        end
        [iv,ix_] = s2mpjlib('ii','46',ix_);
        pb.xnames{iv} = '46';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.0;
        end
        [iv,ix_] = s2mpjlib('ii','46',ix_);
        pb.xnames{iv} = '46';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','47',ix_);
        pb.xnames{iv} = '47';
        ig = ig_('23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','47',ix_);
        pb.xnames{iv} = '47';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -8.88+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -8.88;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.34+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.34;
        end
        [iv,ix_] = s2mpjlib('ii','47',ix_);
        pb.xnames{iv} = '47';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 12.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 12.0;
        end
        [iv,ix_] = s2mpjlib('ii','47',ix_);
        pb.xnames{iv} = '47';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','48',ix_);
        pb.xnames{iv} = '48';
        ig = ig_('19');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','48',ix_);
        pb.xnames{iv} = '48';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.47+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.47;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.46+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.46;
        end
        [iv,ix_] = s2mpjlib('ii','48',ix_);
        pb.xnames{iv} = '48';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','48',ix_);
        pb.xnames{iv} = '48';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .233+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .233;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.358+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.358;
        end
        [iv,ix_] = s2mpjlib('ii','49',ix_);
        pb.xnames{iv} = '49';
        ig = ig_('20');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','49',ix_);
        pb.xnames{iv} = '49';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.07+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.07;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.21+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.21;
        end
        [iv,ix_] = s2mpjlib('ii','49',ix_);
        pb.xnames{iv} = '49';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.5;
        end
        [iv,ix_] = s2mpjlib('ii','49',ix_);
        pb.xnames{iv} = '49';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .205+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .205;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.333+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.333;
        end
        [iv,ix_] = s2mpjlib('ii','50',ix_);
        pb.xnames{iv} = '50';
        ig = ig_('27');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','50',ix_);
        pb.xnames{iv} = '50';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.27+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.27;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.14+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.14;
        end
        [iv,ix_] = s2mpjlib('ii','50',ix_);
        pb.xnames{iv} = '50';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.0;
        end
        [iv,ix_] = s2mpjlib('ii','50',ix_);
        pb.xnames{iv} = '50';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .318+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .318;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.509+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.509;
        end
        [iv,ix_] = s2mpjlib('ii','51',ix_);
        pb.xnames{iv} = '51';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','51',ix_);
        pb.xnames{iv} = '51';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.78+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.78;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.79+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.79;
        end
        [iv,ix_] = s2mpjlib('ii','51',ix_);
        pb.xnames{iv} = '51';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 66.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 66.0;
        end
        [iv,ix_] = s2mpjlib('ii','51',ix_);
        pb.xnames{iv} = '51';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','52',ix_);
        pb.xnames{iv} = '52';
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.426+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.426;
        end
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.204+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.204;
        end
        [iv,ix_] = s2mpjlib('ii','52',ix_);
        pb.xnames{iv} = '52';
        ig = ig_('58');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0924+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0924;
        end
        [iv,ix_] = s2mpjlib('ii','53',ix_);
        pb.xnames{iv} = '53';
        ig = ig_('36');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('56');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.05;
        end
        [iv,ix_] = s2mpjlib('ii','53',ix_);
        pb.xnames{iv} = '53';
        ig = ig_('57');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 9.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 9.05;
        end
        ig = ig_('59');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -9.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -9.5;
        end
        [iv,ix_] = s2mpjlib('ii','53',ix_);
        pb.xnames{iv} = '53';
        ig = ig_('60');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5;
        end
        ig = ig_('61');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .5;
        end
        [iv,ix_] = s2mpjlib('ii','53',ix_);
        pb.xnames{iv} = '53';
        ig = ig_('73');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.36+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.36;
        end
        ig = ig_('74');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.65+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.65;
        end
        [iv,ix_] = s2mpjlib('ii','53',ix_);
        pb.xnames{iv} = '53';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.51+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.51;
        end
        [iv,ix_] = s2mpjlib('ii','54',ix_);
        pb.xnames{iv} = '54';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('26');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','55',ix_);
        pb.xnames{iv} = '55';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('37');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','56',ix_);
        pb.xnames{iv} = '56';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('37');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','57',ix_);
        pb.xnames{iv} = '57';
        ig = ig_('37');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.75+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.75;
        end
        [iv,ix_] = s2mpjlib('ii','58',ix_);
        pb.xnames{iv} = '58';
        ig = ig_('11');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('38');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','58',ix_);
        pb.xnames{iv} = '58';
        ig = ig_('63');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -14.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -14.0;
        end
        ig = ig_('64');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 14.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 14.0;
        end
        [iv,ix_] = s2mpjlib('ii','59',ix_);
        pb.xnames{iv} = '59';
        ig = ig_('12');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('38');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','59',ix_);
        pb.xnames{iv} = '59';
        ig = ig_('63');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.8+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.8;
        end
        ig = ig_('64');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .8+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .8;
        end
        [iv,ix_] = s2mpjlib('ii','60',ix_);
        pb.xnames{iv} = '60';
        ig = ig_('38');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('63');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        [iv,ix_] = s2mpjlib('ii','60',ix_);
        pb.xnames{iv} = '60';
        ig = ig_('64');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.2+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.2;
        end
        [iv,ix_] = s2mpjlib('ii','61',ix_);
        pb.xnames{iv} = '61';
        ig = ig_('4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('34');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','62',ix_);
        pb.xnames{iv} = '62';
        ig = ig_('3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('34');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','63',ix_);
        pb.xnames{iv} = '63';
        ig = ig_('34');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','63',ix_);
        pb.xnames{iv} = '63';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.6+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.6;
        end
        [iv,ix_] = s2mpjlib('ii','64',ix_);
        pb.xnames{iv} = '64';
        ig = ig_('7');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','64',ix_);
        pb.xnames{iv} = '64';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.1;
        end
        [iv,ix_] = s2mpjlib('ii','65',ix_);
        pb.xnames{iv} = '65';
        ig = ig_('8');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','65',ix_);
        pb.xnames{iv} = '65';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 12.63+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 12.63;
        end
        [iv,ix_] = s2mpjlib('ii','66',ix_);
        pb.xnames{iv} = '66';
        ig = ig_('6');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','66',ix_);
        pb.xnames{iv} = '66';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.05;
        end
        ig = ig_('66');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','67',ix_);
        pb.xnames{iv} = '67';
        ig = ig_('5');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','67',ix_);
        pb.xnames{iv} = '67';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 6.9+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 6.9;
        end
        ig = ig_('65');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','68',ix_);
        pb.xnames{iv} = '68';
        ig = ig_('29');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','68',ix_);
        pb.xnames{iv} = '68';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 8.05+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 8.05;
        end
        [iv,ix_] = s2mpjlib('ii','69',ix_);
        pb.xnames{iv} = '69';
        ig = ig_('28');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [iv,ix_] = s2mpjlib('ii','69',ix_);
        pb.xnames{iv} = '69';
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 4.4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 4.4;
        end
        [iv,ix_] = s2mpjlib('ii','70',ix_);
        pb.xnames{iv} = '70';
        ig = ig_('35');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('62');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -10.1+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -10.1;
        end
        [iv,ix_] = s2mpjlib('ii','70',ix_);
        pb.xnames{iv} = '70';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -2.0;
        end
        [iv,ix_] = s2mpjlib('ii','71',ix_);
        pb.xnames{iv} = '71';
        ig = ig_('39');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.325+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.325;
        end
        [iv,ix_] = s2mpjlib('ii','72',ix_);
        pb.xnames{iv} = '72';
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.153+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.153;
        end
        [iv,ix_] = s2mpjlib('ii','73',ix_);
        pb.xnames{iv} = '73';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.316+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.316;
        end
        [iv,ix_] = s2mpjlib('ii','74',ix_);
        pb.xnames{iv} = '74';
        ig = ig_('9');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.814+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.814;
        end
        [iv,ix_] = s2mpjlib('ii','75',ix_);
        pb.xnames{iv} = '75';
        ig = ig_('25');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -3.808+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -3.808;
        end
        [iv,ix_] = s2mpjlib('ii','76',ix_);
        pb.xnames{iv} = '76';
        ig = ig_('24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -4.44+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -4.44;
        end
        [iv,ix_] = s2mpjlib('ii','77',ix_);
        pb.xnames{iv} = '77';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.42+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.42;
        end
        [iv,ix_] = s2mpjlib('ii','77',ix_);
        pb.xnames{iv} = '77';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .04;
        end
        [iv,ix_] = s2mpjlib('ii','78',ix_);
        pb.xnames{iv} = '78';
        ig = ig_('40');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','79',ix_);
        pb.xnames{iv} = '79';
        ig = ig_('10');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5;
        end
        ig = ig_('13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -.5+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -.5;
        end
        [iv,ix_] = s2mpjlib('ii','79',ix_);
        pb.xnames{iv} = '79';
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 3.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 3.0;
        end
        [iv,ix_] = s2mpjlib('ii','80',ix_);
        pb.xnames{iv} = '80';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .4+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .4;
        end
        [iv,ix_] = s2mpjlib('ii','81',ix_);
        pb.xnames{iv} = '81';
        ig = ig_('41');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [iv,ix_] = s2mpjlib('ii','82',ix_);
        pb.xnames{iv} = '82';
        ig = ig_('42');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .0132+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .0132;
        end
        [iv,ix_] = s2mpjlib('ii','83',ix_);
        pb.xnames{iv} = '83';
        ig = ig_('43');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        ig = ig_('C');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = .01+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = .01;
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
        pbm.gconst(ig_('65')) = 23.26;
        pbm.gconst(ig_('66')) = 5.25;
        pbm.gconst(ig_('67')) = 26.32;
        pbm.gconst(ig_('68')) = 21.05;
        pbm.gconst(ig_('69')) = 13.45;
        pbm.gconst(ig_('70')) = 2.58;
        pbm.gconst(ig_('71')) = 10.0;
        pbm.gconst(ig_('72')) = 10.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'D';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.000000000;
        ename = 'E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.109756112;
        ename = 'E3';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.219512224;
        ename = 'E4';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.329268336;
        ename = 'E5';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.439024448;
        ename = 'E6';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.548780441;
        ename = 'E7';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.658536553;
        ename = 'E8';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.768292665;
        ename = 'E9';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.878048778;
        ename = 'E10';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.987804890;
        ename = 'E11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.097560883;
        ename = 'E12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.207317114;
        ename = 'E13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.317073107;
        ename = 'E14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.426829338;
        ename = 'E15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.536585331;
        ename = 'E16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.646341562;
        ename = 'E17';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.756097555;
        ename = 'E18';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.865853548;
        ename = 'E19';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.975609779;
        ename = 'E20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.085365772;
        ename = 'E21';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.195122004;
        ename = 'E22';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.304877996;
        ename = 'E23';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '23';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.414634228;
        ename = 'E24';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '24';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.524390221;
        ename = 'E25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '25';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.634146452;
        ename = 'E26';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '26';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.743902445;
        ename = 'E27';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '27';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.853658438;
        ename = 'E28';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '28';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.963414669;
        ename = 'E29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.073170662;
        ename = 'E30';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '30';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.182926655;
        ename = 'E31';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '31';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.292683125;
        ename = 'E32';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.402439117;
        ename = 'E33';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '33';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.512195110;
        ename = 'E34';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '34';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.621951103;
        ename = 'E35';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '35';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.731707096;
        ename = 'E36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.841463566;
        ename = 'E37';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '37';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 4.951219559;
        ename = 'E38';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '38';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.060975552;
        ename = 'E39';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '39';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.170731544;
        ename = 'E40';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '40';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.280488014;
        ename = 'E41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.390244007;
        ename = 'E42';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '42';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.500000000;
        ename = 'E43';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '43';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.609755993;
        ename = 'E44';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '44';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.719511986;
        ename = 'E45';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '45';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.829268456;
        ename = 'E46';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '46';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 5.939024448;
        ename = 'E47';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '47';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.048780441;
        ename = 'E48';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '48';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.158536434;
        ename = 'E49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.268292904;
        ename = 'E50';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '50';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.378048897;
        ename = 'E51';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '51';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.487804890;
        ename = 'E52';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '52';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.597560883;
        ename = 'E53';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '53';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.707316875;
        ename = 'E54';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '54';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.817073345;
        ename = 'E55';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '55';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 6.926829338;
        ename = 'E56';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '56';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.036585331;
        ename = 'E57';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '57';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.146341324;
        ename = 'E58';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '58';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.256097794;
        ename = 'E59';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '59';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.365853786;
        ename = 'E60';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '60';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.475609779;
        ename = 'E61';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '61';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.585365772;
        ename = 'E62';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '62';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.695121765;
        ename = 'E63';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '63';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.804878235;
        ename = 'E64';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '64';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 7.914634228;
        ename = 'E65';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '65';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.024390221;
        ename = 'E66';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '66';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.134146690;
        ename = 'E67';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '67';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.243902206;
        ename = 'E68';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '68';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.353658676;
        ename = 'E69';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '69';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.463414192;
        ename = 'E70';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '70';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.573170662;
        ename = 'E71';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '71';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.682927132;
        ename = 'E72';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '72';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.792682648;
        ename = 'E73';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '73';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.902439117;
        ename = 'E74';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '74';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.012195587;
        ename = 'E75';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '75';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.121951103;
        ename = 'E76';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '76';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.231707573;
        ename = 'E77';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '77';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.341463089;
        ename = 'E78';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '78';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.451219559;
        ename = 'E79';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '79';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.560976028;
        ename = 'E80';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '80';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.670731544;
        ename = 'E81';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '81';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.780488014;
        ename = 'E82';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '82';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.890243530;
        ename = 'E83';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
        end
        vname = '83';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('D',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 10.000000000;
        v_('1') = 1;
        v_('N') = 83;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('C');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
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
        pb.pbclass = 'C-CQLR2-MN-83-74';
        pb.x0          = zeros(pb.n,1);
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

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = pbm.elpar{iel_}(1)*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0*pbm.elpar{iel_}(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*pbm.elpar{iel_}(1);
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

