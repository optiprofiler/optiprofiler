function varargout = BQPGABIM(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BQPGABIM
%    *********
% 
%    The first 50 variable subproblem from BQPGAUSS.
%    but with variables 1, 15, 42 and 50 fixed at zero
% 
%    SIF input: N. Gould, July 1990.
% 
%    classification = 'C-CQBR2-AN-50-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BQPGABIM';

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
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        irA  = [];
        icA  = [];
        valA = [];
        [ig,ig_] = s2mpjlib('ii','LINGROUP',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        ngrp   = ig_.Count;
        [iv,ix_] = s2mpjlib('ii','1',ix_);
        pb.xnames{iv} = '1';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 5.6987e-02;
        [iv,ix_] = s2mpjlib('ii','2',ix_);
        pb.xnames{iv} = '2';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -6.1847e-03;
        [iv,ix_] = s2mpjlib('ii','3',ix_);
        pb.xnames{iv} = '3';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 5.2516e-03;
        [iv,ix_] = s2mpjlib('ii','4',ix_);
        pb.xnames{iv} = '4';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.1729e-02;
        [iv,ix_] = s2mpjlib('ii','5',ix_);
        pb.xnames{iv} = '5';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 4.9596e-03;
        [iv,ix_] = s2mpjlib('ii','6',ix_);
        pb.xnames{iv} = '6';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.9271e-03;
        [iv,ix_] = s2mpjlib('ii','7',ix_);
        pb.xnames{iv} = '7';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.2185e-02;
        [iv,ix_] = s2mpjlib('ii','8',ix_);
        pb.xnames{iv} = '8';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.3238e-02;
        [iv,ix_] = s2mpjlib('ii','9',ix_);
        pb.xnames{iv} = '9';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.5134e-02;
        [iv,ix_] = s2mpjlib('ii','10',ix_);
        pb.xnames{iv} = '10';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.2247e-02;
        [iv,ix_] = s2mpjlib('ii','11',ix_);
        pb.xnames{iv} = '11';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 2.3741e-02;
        [iv,ix_] = s2mpjlib('ii','12',ix_);
        pb.xnames{iv} = '12';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -9.7666e-02;
        [iv,ix_] = s2mpjlib('ii','13',ix_);
        pb.xnames{iv} = '13';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 9.8702e-02;
        [iv,ix_] = s2mpjlib('ii','14',ix_);
        pb.xnames{iv} = '14';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 7.8901e-04;
        [iv,ix_] = s2mpjlib('ii','15',ix_);
        pb.xnames{iv} = '15';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 5.1663e-04;
        [iv,ix_] = s2mpjlib('ii','16',ix_);
        pb.xnames{iv} = '16';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.7477e-04;
        [iv,ix_] = s2mpjlib('ii','17',ix_);
        pb.xnames{iv} = '17';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.1795e-03;
        [iv,ix_] = s2mpjlib('ii','18',ix_);
        pb.xnames{iv} = '18';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.7351e-02;
        [iv,ix_] = s2mpjlib('ii','19',ix_);
        pb.xnames{iv} = '19';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.3439e-03;
        [iv,ix_] = s2mpjlib('ii','20',ix_);
        pb.xnames{iv} = '20';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -5.6977e-02;
        [iv,ix_] = s2mpjlib('ii','21',ix_);
        pb.xnames{iv} = '21';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.0040e-02;
        [iv,ix_] = s2mpjlib('ii','22',ix_);
        pb.xnames{iv} = '22';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -8.3380e-02;
        [iv,ix_] = s2mpjlib('ii','23',ix_);
        pb.xnames{iv} = '23';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -3.7526e-03;
        [iv,ix_] = s2mpjlib('ii','24',ix_);
        pb.xnames{iv} = '24';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -9.4555e-04;
        [iv,ix_] = s2mpjlib('ii','25',ix_);
        pb.xnames{iv} = '25';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.9258e-03;
        [iv,ix_] = s2mpjlib('ii','26',ix_);
        pb.xnames{iv} = '26';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.3959e-03;
        [iv,ix_] = s2mpjlib('ii','27',ix_);
        pb.xnames{iv} = '27';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.3749e-03;
        [iv,ix_] = s2mpjlib('ii','28',ix_);
        pb.xnames{iv} = '28';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.3677e-03;
        [iv,ix_] = s2mpjlib('ii','29',ix_);
        pb.xnames{iv} = '29';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -2.7985e-02;
        [iv,ix_] = s2mpjlib('ii','30',ix_);
        pb.xnames{iv} = '30';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.8839e-03;
        [iv,ix_] = s2mpjlib('ii','31',ix_);
        pb.xnames{iv} = '31';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.2340e-03;
        [iv,ix_] = s2mpjlib('ii','32',ix_);
        pb.xnames{iv} = '32';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -6.8139e-04;
        [iv,ix_] = s2mpjlib('ii','33',ix_);
        pb.xnames{iv} = '33';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -3.5838e-02;
        [iv,ix_] = s2mpjlib('ii','34',ix_);
        pb.xnames{iv} = '34';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -3.4857e-02;
        [iv,ix_] = s2mpjlib('ii','35',ix_);
        pb.xnames{iv} = '35';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 2.8724e-03;
        [iv,ix_] = s2mpjlib('ii','36',ix_);
        pb.xnames{iv} = '36';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.6625e-02;
        [iv,ix_] = s2mpjlib('ii','37',ix_);
        pb.xnames{iv} = '37';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 1.3571e-02;
        [iv,ix_] = s2mpjlib('ii','38',ix_);
        pb.xnames{iv} = '38';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -7.2447e-03;
        [iv,ix_] = s2mpjlib('ii','39',ix_);
        pb.xnames{iv} = '39';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.6034e-04;
        [iv,ix_] = s2mpjlib('ii','40',ix_);
        pb.xnames{iv} = '40';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.6225e-02;
        [iv,ix_] = s2mpjlib('ii','41',ix_);
        pb.xnames{iv} = '41';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 2.2034e-05;
        [iv,ix_] = s2mpjlib('ii','42',ix_);
        pb.xnames{iv} = '42';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 5.8844e-02;
        [iv,ix_] = s2mpjlib('ii','43',ix_);
        pb.xnames{iv} = '43';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 3.0725e-03;
        [iv,ix_] = s2mpjlib('ii','44',ix_);
        pb.xnames{iv} = '44';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 2.8227e-03;
        [iv,ix_] = s2mpjlib('ii','45',ix_);
        pb.xnames{iv} = '45';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -2.0681e-02;
        [iv,ix_] = s2mpjlib('ii','46',ix_);
        pb.xnames{iv} = '46';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -5.4952e-03;
        [iv,ix_] = s2mpjlib('ii','47',ix_);
        pb.xnames{iv} = '47';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 6.2552e-04;
        [iv,ix_] = s2mpjlib('ii','48',ix_);
        pb.xnames{iv} = '48';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = 3.3782e-02;
        [iv,ix_] = s2mpjlib('ii','49',ix_);
        pb.xnames{iv} = '49';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -4.8584e-03;
        [iv,ix_] = s2mpjlib('ii','50',ix_);
        pb.xnames{iv} = '50';
        icA(end+1) = iv;
        irA(end+1) = ig_('LINGROUP');
        valA(end+1) = -1.4371e-03;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -0.1*ones(pb.n,1);
        pb.xupper = 0.1*ones(pb.n,1);
        pb.xlower(ix_('1'),1) = 0.0;
        pb.xupper(ix_('1'),1) = 0.0;
        pb.xlower(ix_('2'),1) = -3.9206e-03;
        pb.xupper(ix_('3')) = 9.9999e-02;
        pb.xlower(ix_('4'),1) = -1.0001e-01;
        pb.xupper(ix_('4')) = 9.9990e-02;
        pb.xupper(ix_('5')) = 9.9997e-02;
        pb.xlower(ix_('6'),1) = -9.9994e-02;
        pb.xupper(ix_('6')) = 6.1561e-06;
        pb.xlower(ix_('7'),1) = -3.9119e-03;
        pb.xupper(ix_('7')) = 9.9986e-02;
        pb.xlower(ix_('8'),1) = -1.0001e-01;
        pb.xupper(ix_('8')) = 2.5683e-02;
        pb.xlower(ix_('9'),1) = -9.9987e-02;
        pb.xupper(ix_('9')) = 1.0001e-01;
        pb.xlower(ix_('10'),1) = -9.9988e-02;
        pb.xupper(ix_('10')) = 1.0001e-01;
        pb.xlower(ix_('11'),1) = -1.0001e-01;
        pb.xupper(ix_('11')) = 2.8998e-03;
        pb.xlower(ix_('12'),1) = -9.9952e-02;
        pb.xupper(ix_('12')) = 4.7652e-05;
        pb.xlower(ix_('13'),1) = -4.5551e-05;
        pb.xupper(ix_('13')) = 9.9954e-02;
        pb.xlower(ix_('14'),1) = -9.9999e-02;
        pb.xlower(ix_('15'),1) = 0.0;
        pb.xupper(ix_('15'),1) = 0.0;
        pb.xlower(ix_('16'),1) = -7.2801e-02;
        pb.xlower(ix_('18'),1) = -9.9992e-02;
        pb.xupper(ix_('18')) = 8.3681e-06;
        pb.xlower(ix_('20'),1) = -9.9956e-02;
        pb.xupper(ix_('20')) = 4.3809e-05;
        pb.xlower(ix_('22'),1) = -9.9961e-02;
        pb.xupper(ix_('22')) = 3.9248e-05;
        pb.xlower(ix_('25'),1) = -4.1110e-03;
        pb.xlower(ix_('29'),1) = -9.6988e-02;
        pb.xupper(ix_('29')) = 1.0002e-01;
        pb.xlower(ix_('32'),1) = -5.8439e-02;
        pb.xlower(ix_('33'),1) = -4.5616e-06;
        pb.xupper(ix_('33')) = 9.9995e-02;
        pb.xlower(ix_('34'),1) = -9.9999e-02;
        pb.xupper(ix_('34')) = 7.3117e-07;
        pb.xlower(ix_('35'),1) = -9.9991e-02;
        pb.xupper(ix_('35')) = 9.3168e-06;
        pb.xlower(ix_('36'),1) = -9.9977e-02;
        pb.xupper(ix_('36')) = 1.0002e-01;
        pb.xlower(ix_('37'),1) = -9.9984e-02;
        pb.xupper(ix_('37')) = 1.5812e-05;
        pb.xlower(ix_('39'),1) = -3.9611e-06;
        pb.xupper(ix_('39')) = 9.9996e-02;
        pb.xlower(ix_('40'),1) = -8.8262e-06;
        pb.xupper(ix_('40')) = 9.9991e-02;
        pb.xlower(ix_('41'),1) = -1.0001e-01;
        pb.xupper(ix_('41')) = 9.9986e-02;
        pb.xlower(ix_('42'),1) = 0.0;
        pb.xupper(ix_('42'),1) = 0.0;
        pb.xlower(ix_('43'),1) = -1.9873e-06;
        pb.xupper(ix_('43')) = 9.9998e-02;
        pb.xlower(ix_('45'),1) = -9.9993e-02;
        pb.xupper(ix_('45')) = 7.4220e-06;
        pb.xlower(ix_('46'),1) = -9.9999e-02;
        pb.xupper(ix_('46')) = 8.2308e-07;
        pb.xlower(ix_('47'),1) = -3.0424e-06;
        pb.xupper(ix_('47')) = 9.9997e-02;
        pb.xlower(ix_('48'),1) = -9.9985e-02;
        pb.xupper(ix_('48')) = 1.5119e-05;
        pb.xlower(ix_('49'),1) = -1.0004e-01;
        pb.xupper(ix_('49')) = 2.4305e-02;
        pb.xlower(ix_('50'),1) = 0.0;
        pb.xupper(ix_('50'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOFFDIAG',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eDIAG',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'D   1   1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  11  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  12  12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '12';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  20  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  21';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  21';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  21  21';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '21';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  29  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  36  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  37';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '37';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  36  37';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '37';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  37  37';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '37';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  41  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  42';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '42';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  42';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '42';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  42  42';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '42';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   1  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  49  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   2   2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  13  13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  22';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  22';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  22  22';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '22';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  30';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '30';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  30';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '30';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  30  30';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '30';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  38';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '38';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  36  38';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '38';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  38  38';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '38';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   2  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   3   3';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  23';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '23';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  23';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '23';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  23  23';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '23';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  43';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '43';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  43';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '43';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  43  43';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '43';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   3  50';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '50';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  49  50';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '50';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  50  50';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '50';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   4   4';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  14  14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  31';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '31';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  31';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '31';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  31  31';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '31';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  44';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '44';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  44';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '44';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  44  44';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '44';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   4  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   5   5';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  15  15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '15';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  24';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '24';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  24';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '24';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  24  24';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '24';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  32';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  32';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  32  32';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '32';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  39';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '39';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  36  39';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '39';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  39  39';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '39';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  45';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '45';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  45';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '45';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  45  45';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '45';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   5  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   6   6';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  16  16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '16';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '25';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '25';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  25  25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '25';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  33';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '33';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  33';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '33';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  33  33';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '33';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  46';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '46';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  46';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '46';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  46  46';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '46';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   6  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   7   7';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  17';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  17';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  17  17';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '17';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  34';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '34';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  34';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '34';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  34  34';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '34';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   7  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   8   8';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  40';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '40';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  36  40';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '40';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  40  40';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '40';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   8  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D   9   9';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  26';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '26';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  26';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '26';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  26  26';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '26';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  35';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '35';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  29  35';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '35';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  35  35';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '35';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O   9  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  10  10';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  29';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '29';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  36';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '36';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  41';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  47';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '47';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  47';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '47';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  47  47';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '47';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  10  49';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '49';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  18';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  18  18';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '18';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  27';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '27';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  27  27';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '27';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  11  19';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  19  19';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '19';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  20  28';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '20';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '28';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  28  28';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '28';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'O  41  48';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eOFFDIAG';
            ielftype(ie) = iet_('eOFFDIAG');
        end
        vname = '41';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = '48';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'D  48  48';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eDIAG';
        ielftype(ie) = iet_('eDIAG');
        vname = '48';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.1,0.1,[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('LINGROUP');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   1   1');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  11');
        pbm.grelw{ig}(posel) = -9.9819e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  11  11');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  12');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  12');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  12  12');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  20');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  20  20');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  21');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  21');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  21  21');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  29');
        pbm.grelw{ig}(posel) = 9.0362e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  29  29');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  36');
        pbm.grelw{ig}(posel) = 6.5103e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  36  36');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  37');
        pbm.grelw{ig}(posel) = 6.5140e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  36  37');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  37  37');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  41');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  41  41');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  42');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  42');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  42  42');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   1  49');
        pbm.grelw{ig}(posel) = -9.7537e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  49  49');
        pbm.grelw{ig}(posel) = 7.8331e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   2   2');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  11');
        pbm.grelw{ig}(posel) = -9.9213e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  13');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  13');
        pbm.grelw{ig}(posel) = 9.9608e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  13  13');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  20');
        pbm.grelw{ig}(posel) = -9.9698e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  22');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  22');
        pbm.grelw{ig}(posel) = 9.9608e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  22  22');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  29');
        pbm.grelw{ig}(posel) = 8.9945e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  30');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  30');
        pbm.grelw{ig}(posel) = 9.9608e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  30  30');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  36');
        pbm.grelw{ig}(posel) = 6.4885e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  38');
        pbm.grelw{ig}(posel) = 6.5140e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  36  38');
        pbm.grelw{ig}(posel) = 9.9608e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  38  38');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  41');
        pbm.grelw{ig}(posel) = 7.5197e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   2  49');
        pbm.grelw{ig}(posel) = -9.7167e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   3   3');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  11');
        pbm.grelw{ig}(posel) = 8.1209e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  20');
        pbm.grelw{ig}(posel) = 8.1463e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  23');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  23');
        pbm.grelw{ig}(posel) = -8.1463e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  23  23');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  29');
        pbm.grelw{ig}(posel) = -7.3536e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  36');
        pbm.grelw{ig}(posel) = -5.3119e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  41');
        pbm.grelw{ig}(posel) = -6.1506e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  43');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  43');
        pbm.grelw{ig}(posel) = -8.1463e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  43  43');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  49');
        pbm.grelw{ig}(posel) = 7.9480e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   3  50');
        pbm.grelw{ig}(posel) = -9.7566e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  49  50');
        pbm.grelw{ig}(posel) = -8.1463e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  50  50');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   4   4');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  11');
        pbm.grelw{ig}(posel) = 2.8141e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  14');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  14');
        pbm.grelw{ig}(posel) = -2.8225e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  14  14');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  20');
        pbm.grelw{ig}(posel) = 2.8228e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  29');
        pbm.grelw{ig}(posel) = -2.5487e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  31');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  31');
        pbm.grelw{ig}(posel) = -2.8225e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  31  31');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  36');
        pbm.grelw{ig}(posel) = -1.8370e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  41');
        pbm.grelw{ig}(posel) = -2.1312e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  44');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  44');
        pbm.grelw{ig}(posel) = -2.8225e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  44  44');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   4  49');
        pbm.grelw{ig}(posel) = 2.7539e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   5   5');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  11');
        pbm.grelw{ig}(posel) = 2.6350e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  15');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  15');
        pbm.grelw{ig}(posel) = -2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  15  15');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  20');
        pbm.grelw{ig}(posel) = 2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  24');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  24');
        pbm.grelw{ig}(posel) = -2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  24  24');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  29');
        pbm.grelw{ig}(posel) = -2.3863e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  32');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  32');
        pbm.grelw{ig}(posel) = -2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  32  32');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  36');
        pbm.grelw{ig}(posel) = -1.7205e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  39');
        pbm.grelw{ig}(posel) = 6.5140e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  36  39');
        pbm.grelw{ig}(posel) = -2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  39  39');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  41');
        pbm.grelw{ig}(posel) = -1.9971e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  45');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  45');
        pbm.grelw{ig}(posel) = -2.6427e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  45  45');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   5  49');
        pbm.grelw{ig}(posel) = 2.5757e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   6   6');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  11');
        pbm.grelw{ig}(posel) = 9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  16');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  16');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  16  16');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  20');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  25');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  25');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  25  25');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  29');
        pbm.grelw{ig}(posel) = -9.0289e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  33');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  33');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  33  33');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  36');
        pbm.grelw{ig}(posel) = -6.5144e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  41');
        pbm.grelw{ig}(posel) = -7.5509e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  46');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  46');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  46  46');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   6  49');
        pbm.grelw{ig}(posel) = 9.7565e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   7   7');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  11');
        pbm.grelw{ig}(posel) = -9.9320e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  17');
        pbm.grelw{ig}(posel) = -9.9709e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  17');
        pbm.grelw{ig}(posel) = 9.9610e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  17  17');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  20');
        pbm.grelw{ig}(posel) = -9.9631e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  29');
        pbm.grelw{ig}(posel) = 8.9946e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  34');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  34');
        pbm.grelw{ig}(posel) = 9.9610e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  34  34');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  36');
        pbm.grelw{ig}(posel) = 6.4890e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  41');
        pbm.grelw{ig}(posel) = 7.5199e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   7  49');
        pbm.grelw{ig}(posel) = -9.7188e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   8   8');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  11');
        pbm.grelw{ig}(posel) = 9.7157e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  20');
        pbm.grelw{ig}(posel) = 9.7417e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  29');
        pbm.grelw{ig}(posel) = -8.7973e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  36');
        pbm.grelw{ig}(posel) = -6.3446e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  40');
        pbm.grelw{ig}(posel) = 6.5140e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  36  40');
        pbm.grelw{ig}(posel) = -9.7431e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  40  40');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  41');
        pbm.grelw{ig}(posel) = -7.3586e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   8  49');
        pbm.grelw{ig}(posel) = 9.5052e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D   9   9');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  11');
        pbm.grelw{ig}(posel) = -2.9055;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  20');
        pbm.grelw{ig}(posel) = -2.9605;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  26');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  26');
        pbm.grelw{ig}(posel) = 2.9604;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  26  26');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  29');
        pbm.grelw{ig}(posel) = 2.6517;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  35');
        pbm.grelw{ig}(posel) = 9.0300e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  29  35');
        pbm.grelw{ig}(posel) = 2.9604;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  35  35');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  36');
        pbm.grelw{ig}(posel) = 1.9168;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  41');
        pbm.grelw{ig}(posel) = 2.2464;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O   9  49');
        pbm.grelw{ig}(posel) = -2.9243;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  10  10');
        pbm.grelw{ig}(posel) = 1.0624e+03;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  11');
        pbm.grelw{ig}(posel) = 2.9135e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  20');
        pbm.grelw{ig}(posel) = 2.9241e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  29');
        pbm.grelw{ig}(posel) = -2.6379e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  36');
        pbm.grelw{ig}(posel) = -1.9046e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  41');
        pbm.grelw{ig}(posel) = -2.2065e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  10  47');
        pbm.grelw{ig}(posel) = 7.5507e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  47');
        pbm.grelw{ig}(posel) = -2.9232e+01;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  47  47');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  18');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  18  18');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  27');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  27  27');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  11  19');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  19  19');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  20  28');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  28  28');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('O  41  48');
        pbm.grelw{ig}(posel) = -1.0000e+02;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('D  48  48');
        pbm.grelw{ig}(posel) = 1.0000e+02;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO BQPGABIM             -3.790343D-5
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CQBR2-AN-50-0';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

    case 'eDIAG'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 0.5*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0;
                varargout{3} = H_;
            end
        end

    case 'eOFFDIAG'

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

