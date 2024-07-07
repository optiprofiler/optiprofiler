function varargout = SANTALS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SANTALS
%    --------
% 
%    The Santa problem as suggested in a Christmas competition
%    by Jens Jensen (Scientific Computing, STFC). To quote Jens,
% 
%    Santa and His Elves
% 
%    SCD Christmas programming challenge 2016
% 
%    Christmas has come to the Santa Claus Department – or rather, the SCD
%    is coming to Christmas. Santa is flying around the world, presently
%    presenting presents. Ho, ho, ho! No striking air crew on Santa’s sleigh!
%    No airport strikes on the North Pole.
% 
%    For the purpose of this exercise, the Earth is round as a perfect ball,
%    with radius precisely 6,371,000 metres. However, everything is at the
%    same longitude and latitude as the “real” Earth. So for example, the
%    Greenwich observatory is at 51°28'40"N 0°00'04"W both on the “real”
%    Earth and on Santa’s Earth. (Also ignore rotation of the Earth and
%    anything practical like that.)
% 
%    Santa sets off from the North Pole along 2°6'57.6" E bearing south
%    (obviously), and also bearing presents (obviously). Whenever Santa
%    leaves a location, he leaves an elf behind, in order to help unwrapping
%    presents; the elf does this and then flies out independently to meet up
%    with Santa at the next location - this means that Santa only needs two
%    Elves. Here’s how:
% 
%    1. Santa leaves the North Pole, setting out for location A. Elf 1 is
%    left behind (in this particular case, not to unwrap presents, but to
%    turn the lights off, and ensure the oven is off – it's elf'n'safety,
%    you know.)
% 
%    2. Santa arrives in location A and hands them their present. Now Elf 2
%    is with Santa; Elf 1 is still at the NP.
% 
%    3. Santa leaves location A, leaving behind Elf 2. Santa flies on to
%    location B; Elf 1, who remained at the North Pole, also flies to B and
%    meets Santa there; Elf 2 is left behind at A.
% 
%    4. Santa arrives at location B along with Elf 1, and hands out
%    presents. Santa then leaves location B to fly to C, leaving behind Elf 1
%    at location B. Meanwhile Elf 2, having finished helping at location A,
%    leaves location A to fly on to C, to meet Santa there.
% 
%    5. Santa arrives from B at location C; Elf 2 also arrives into C from
%    location A. Elf 1 remains at B until Santa flies onward to location D.
% 
%    6. At the last hop, Santa needs a rest and flies to 31°46'42.4" S
%    144°46'12.9" W.  The Elves also fly to this location - maps show no land
%    here but it is hidden. Either that or we got the coordinates wrong.
%    In either case Santa and elves do fly to this location.
% 
%    The following table shows the distance of Santa's hops, as well as those
%    of the elves, with the distance given in metres:
% 
%    Who     Hop  Distance travelled
%    Santa   1    5405238
%            2    623852
%            3    1005461
%            4    7470967
%            5    3632559
%            6    10206818
%            7    7967212
%            8    5896361
%            9    8337266
%            10   13019505
%            11   8690818
%            12   8971302
%    Elf1    1    4866724
%            2    6833740
%            3    13489586
%            4    9195575
%            5    9704793
%            6    12498127
%    Elf2    1    1375828
%            2    4917407
%            3    10617953
%            4    10996150
%            5    7901038
%            6    8971302
% 
%    What is Santa’s route?  What sort of presents is he carrying?
% 
%    Bonus question: did you really need to know the starting direction?
% 
%    Added by Nick: the problem has many local minimizers, but it is only
%    the global minimizer that is of interest.
% 
%    least-squares version
% 
%    SIF input: Nick Gould, Dec 2016.
% 
%    classification = 'SBR2-AN-21-0'
% 
%    Number of stops on Santa's path (path goes from index 0 to 12)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SANTALS';

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
        v_('S') = 12;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('180.0') = 180.0;
        v_('S-1') = -1+v_('S');
        v_('RS') = v_('S');
        v_('PI/4') = atan(1.0);
        v_('PI') = 4.0*v_('PI/4');
        v_('PI/180') = v_('PI')/v_('180.0');
        v_('PHI0') = 90.0;
        v_('LAM0') = 0.0;
        v_('PHI12') = -31.77844444;
        v_('LAM12') = -144.77025;
        v_('LAM1') = 2.116;
        v_('PHI0') = v_('PHI0')*v_('PI/180');
        v_('LAM0') = v_('LAM0')*v_('PI/180');
        v_('PHI12') = v_('PHI12')*v_('PI/180');
        v_('LAM12') = v_('LAM12')*v_('PI/180');
        v_('LAM1') = v_('LAM1')*v_('PI/180');
        v_('DPHI') = v_('PHI12')-v_('PHI0');
        v_('DLAM') = v_('LAM12')-v_('LAM0');
        v_('DPHI/S') = v_('DPHI')/v_('RS');
        v_('DLAM/S') = v_('DLAM')/v_('RS');
        v_('RADIUS') = 6371000.0;
        v_('D0,1') = 5405238.0;
        v_('D0,2') = 4866724.0;
        v_('D1,2') = 623852.0;
        v_('D1,3') = 1375828.0;
        v_('D2,3') = 1005461.0;
        v_('D2,4') = 6833740.0;
        v_('D3,4') = 7470967.0;
        v_('D3,5') = 4917407.0;
        v_('D4,5') = 3632559.0;
        v_('D4,6') = 13489586.0;
        v_('D5,6') = 10206818.0;
        v_('D5,7') = 10617953.0;
        v_('D6,7') = 7967212.0;
        v_('D6,8') = 9195575.0;
        v_('D7,8') = 5896361.0;
        v_('D7,9') = 10996150.0;
        v_('D8,9') = 8337266.0;
        v_('D8,10') = 9704793.0;
        v_('D9,10') = 13019505.0;
        v_('D9,11') = 7901038.0;
        v_('D10,11') = 8690818.0;
        v_('D10,12') = 12498127.0;
        v_('D11,12') = 8971302.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','PHI1',ix_);
        pb.xnames{iv} = 'PHI1';
        for I=v_('2'):v_('S-1')
            [iv,ix_] = s2mpjlib('ii',['PHI',int2str(I)],ix_);
            pb.xnames{iv} = ['PHI',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['LAM',int2str(I)],ix_);
            pb.xnames{iv} = ['LAM',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','R0,1',ig_);
        gtype{ig} = '<>';
        for I=v_('2'):v_('S')
            v_('I1') = -1+I;
            v_('I2') = -2+I;
            [ig,ig_] = s2mpjlib('ii',['R',int2str(round(v_('I2'))),',',int2str(I)],ig_);
            gtype{ig} = '<>';
            [ig,ig_] = s2mpjlib('ii',['R',int2str(round(v_('I1'))),',',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        v_('D/RAD') = v_('D0,1')/v_('RADIUS');
        v_('CD/RAD') = cos(v_('D/RAD'));
        pbm.gconst(ig_('R0,1')) = v_('CD/RAD');
        for I=v_('2'):v_('S')
            v_('I2') = -2+I;
            v_('D/RAD') =...
                  v_(['D',int2str(round(v_('I2'))),',',int2str(I)])/v_('RADIUS');
            v_('CD/RAD') = cos(v_('D/RAD'));
            pbm.gconst(ig_(['R',int2str(round(v_('I2'))),',',int2str(I)])) =...
                  v_('CD/RAD');
            v_('I1') = -1+I;
            v_('D/RAD') =...
                  v_(['D',int2str(round(v_('I1'))),',',int2str(I)])/v_('RADIUS');
            v_('CD/RAD') = cos(v_('D/RAD'));
            pbm.gconst(ig_(['R',int2str(round(v_('I1'))),',',int2str(I)])) =...
                  v_('CD/RAD');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1000.0*ones(pb.n,1);
        pb.xupper = 1000.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('PHI1'),1) = 0.7223835215;
        pb.x0(ix_('PHI2'),1) = 0.8069093428;
        pb.x0(ix_('LAM2'),1) = -0.031657133;
        pb.x0(ix_('PHI3'),1) = 0.9310164154;
        pb.x0(ix_('LAM3'),1) = 0.1199353230;
        pb.x0(ix_('PHI4'),1) = 6.6067392710;
        pb.x0(ix_('LAM4'),1) = -1.214314477;
        pb.x0(ix_('PHI5'),1) = -3.530946794;
        pb.x0(ix_('LAM5'),1) = 2.5329493980;
        pb.x0(ix_('PHI6'),1) = -9.798251905;
        pb.x0(ix_('LAM6'),1) = 4.3021328700;
        pb.x0(ix_('PHI7'),1) = 14.632267534;
        pb.x0(ix_('LAM7'),1) = -12.96253311;
        pb.x0(ix_('PHI8'),1) = 2.0349445303;
        pb.x0(ix_('LAM8'),1) = -4.050000443;
        pb.x0(ix_('PHI9'),1) = -28.45607804;
        pb.x0(ix_('LAM9'),1) = 22.430117198;
        pb.x0(ix_('PHI10'),1) = 16.034035489;
        pb.x0(ix_('LAM10'),1) = -17.28050167;
        pb.x0(ix_('PHI11'),1) = 0.8717052037;
        pb.x0(ix_('LAM11'),1) = -0.833052840;
        v_('PHIS') = v_('DPHI/S');
        v_('START') = v_('PHI0')+v_('PHIS');
        for I=v_('2'):v_('S-1')
            v_('RI') = I;
            v_('PHIS') = v_('DPHI/S')*v_('RI');
            v_('START') = v_('PHI0')+v_('PHIS');
            v_('LAMS') = v_('DLAM/S')*v_('RI');
            v_('START') = v_('LAM0')+v_('LAMS');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE',iet_);
        elftv{it}{1} = 'PHI1';
        elftv{it}{2} = 'PHI2';
        elftv{it}{3} = 'LAM1';
        elftv{it}{4} = 'LAM2';
        [it,iet_] = s2mpjlib( 'ii', 'eE3',iet_);
        elftv{it}{1} = 'PHI1';
        elftv{it}{2} = 'PHI2';
        elftv{it}{3} = 'LAM1';
        elftp{it}{1} = 'LAMF';
        [it,iet_] = s2mpjlib( 'ii', 'eE2',iet_);
        elftv{it}{1} = 'PHI1';
        elftv{it}{2} = 'LAM1';
        elftp{it}{1} = 'PHIF';
        elftp{it}{2} = 'LAMF';
        [it,iet_] = s2mpjlib( 'ii', 'eE1',iet_);
        elftv{it}{1} = 'PHI1';
        elftp{it}{1} = 'PHIF';
        elftp{it}{2} = 'LAMF';
        elftp{it}{3} = 'LAMS';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E0,1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE1';
        ielftype(ie) = iet_('eE1');
        vname = 'PHI1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PHIF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('PHI0');
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM0');
        [~,posep] = ismember('LAMS',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM1');
        ename = 'E0,2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE2';
        ielftype(ie) = iet_('eE2');
        vname = 'PHI2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PHIF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('PHI0');
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM0');
        ename = 'E1,2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE3';
        ielftype(ie) = iet_('eE3');
        vname = 'PHI1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PHI2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM1');
        ename = 'E1,3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE3';
        ielftype(ie) = iet_('eE3');
        vname = 'PHI1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PHI3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM1');
        ename = 'E2,3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'PHI2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'PHI3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        for I=v_('4'):v_('S-1')
            v_('I2') = -2+I;
            ename = ['E',int2str(round(v_('I2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE';
            ielftype(ie) = iet_('eE');
            ename = ['E',int2str(round(v_('I2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['PHI',int2str(round(v_('I2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['PHI',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('PHI2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['LAM',int2str(round(v_('I2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I2'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['LAM',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('LAM2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I1') = -1+I;
            ename = ['E',int2str(round(v_('I1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE';
            ielftype(ie) = iet_('eE');
            ename = ['E',int2str(round(v_('I1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['PHI',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['PHI',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('PHI2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['LAM',int2str(round(v_('I1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('I1'))),',',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['LAM',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
            posev = find(strcmp('LAM2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        ename = 'E10,12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE2';
        ielftype(ie) = iet_('eE2');
        vname = 'PHI10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PHIF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('PHI12');
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM12');
        ename = 'E11,12';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE2';
        ielftype(ie) = iet_('eE2');
        vname = 'PHI11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('PHI1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'LAM11';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1000.0,1000.0,[]);
        posev = find(strcmp('LAM1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('PHIF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('PHI12');
        [~,posep] = ismember('LAMF',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('LAM12');
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        ig = ig_('R0,1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E0,1');
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('2'):v_('S')
            v_('I2') = -2+I;
            ig = ig_(['R',int2str(round(v_('I2'))),',',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('I2'))),',',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            v_('I1') = -1+I;
            ig = ig_(['R',int2str(round(v_('I1'))),',',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('I1'))),',',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SANTA               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-AN-21-0';
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

    case 'eE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S1 = sin(EV_(1));
        S2 = sin(EV_(2));
        C1 = cos(EV_(1));
        C2 = cos(EV_(2));
        C = cos(EV_(3)-EV_(4));
        S = sin(EV_(3)-EV_(4));
        C1C2S = C1*C2*S;
        C1C2C = C1*C2*C;
        C1S2S = C1*S2*S;
        S1C2S = S1*C2*S;
        varargout{1} = S1*S2+C1*C2*C;
        if(nargout>1)
            g_(1,1) = C1*S2-S1*C2*C;
            g_(2,1) = S1*C2-C1*S2*C;
            g_(3,1) = -C1C2S;
            g_(4,1) = C1C2S;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) = -S1*S2-C1*C2*C;
                H_(2,1) = C1*C2+S1*S2*C;
                H_(1,2) = H_(2,1);
                H_(2,2) = -S1*S2-C1*C2*C;
                H_(3,1) = S1C2S;
                H_(1,3) = H_(3,1);
                H_(3,2) = C1S2S;
                H_(2,3) = H_(3,2);
                H_(3,3) = -C1C2C;
                H_(1,4) = -S1C2S;
                H_(4,1) = H_(1,4);
                H_(2,4) = -C1S2S;
                H_(4,2) = H_(2,4);
                H_(3,4) = C1C2C;
                H_(4,3) = H_(3,4);
                H_(4,4) = -C1C2C;
                varargout{3} = H_;
            end
        end

    case 'eE3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S1 = sin(EV_(1));
        S2 = sin(EV_(2));
        C1 = cos(EV_(1));
        C2 = cos(EV_(2));
        C = cos(EV_(3)-pbm.elpar{iel_}(1));
        S = sin(EV_(3)-pbm.elpar{iel_}(1));
        varargout{1} = S1*S2+C1*C2*C;
        if(nargout>1)
            g_(1,1) = C1*S2-S1*C2*C;
            g_(2,1) = S1*C2-C1*S2*C;
            g_(3,1) = -C1*C2*S;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = -S1*S2-C1*C2*C;
                H_(2,1) = C1*C2+S1*S2*C;
                H_(1,2) = H_(2,1);
                H_(2,2) = -S1*S2-C1*C2*C;
                H_(3,1) = S1*C2*S;
                H_(1,3) = H_(3,1);
                H_(3,2) = C1*S2*S;
                H_(2,3) = H_(3,2);
                H_(3,3) = -C1*C2*C;
                varargout{3} = H_;
            end
        end

    case 'eE2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S1 = sin(EV_(1));
        SF = sin(pbm.elpar{iel_}(1));
        C1 = cos(EV_(1));
        CF = cos(pbm.elpar{iel_}(1));
        C = cos(EV_(2)-pbm.elpar{iel_}(2));
        S = sin(EV_(2)-pbm.elpar{iel_}(2));
        varargout{1} = S1*SF+C1*CF*C;
        if(nargout>1)
            g_(1,1) = C1*SF-S1*CF*C;
            g_(2,1) = -C1*CF*S;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = -S1*SF-C1*CF*C;
                H_(2,1) = S1*CF*S;
                H_(1,2) = H_(2,1);
                H_(2,2) = -C1*CF*C;
                varargout{3} = H_;
            end
        end

    case 'eE1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        S1 = sin(EV_(1));
        SF = sin(pbm.elpar{iel_}(1));
        C1 = cos(EV_(1));
        CF = cos(pbm.elpar{iel_}(1));
        C = cos(pbm.elpar{iel_}(3)-pbm.elpar{iel_}(2));
        S = sin(pbm.elpar{iel_}(3)-pbm.elpar{iel_}(2));
        varargout{1} = S1*SF+C1*CF*C;
        if(nargout>1)
            g_(1,1) = C1*SF-S1*CF*C;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -S1*SF-C1*CF*C;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
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

