function varargout = BAmL1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BAmL1
%    *********
% 
%    Bundle Adjustment problem from reconstructive geometry in which
%    a collection of photographs is used to determine the position of
%    a set of observed points. Each observed point is seen via its
%    two-dimensional projections on a subset of the photographs. The
%    solution is found by solvng a large nonlinear least-squares problem.
%    This variant is given as an inconsistent set of nonlinear equations.
% 
%    Source: data from the Bundle Adjustment in the Large
%    project, http://grail.cs.washington.edu/projects/bal/
% 
%    Ladybug datasets (single image extracted)
% 
%    SIF input: Nick Gould, November 2016
% 
%    classification = 'NOR2-MN-57-12'
% 
%    Number of images
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BAmL1';

switch(action)

    case 'setup'

    pb.name      = 'BAmL1';
    pb.sifpbname = 'BAmL1';
    pbm.name     = 'BAmL1';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('nuimages') = 49;
        v_('nupoints') = 1;
        v_('nuobservs') = 6;
        v_('1') = 1;
        v_('O1') = 1.0;
        v_('O2') = 2.0;
        v_('O3') = 4.0;
        v_('O4') = 27.0;
        v_('O5') = 30.0;
        v_('O6') = 37.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('nupoints')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2xlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        for J=v_('1'):v_('nuobservs')
            v_('RI') = v_(['O',int2str(J)]);
            v_('I') = fix(v_('RI'));
            [iv,ix_] = s2xlib('ii',['RX',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['RX',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['RY',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['RY',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['RZ',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['RZ',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['TX',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['TX',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['TY',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['TY',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['TZ',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['TZ',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['KA',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['KA',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['KB',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['KB',int2str(round(v_('I')))];
            [iv,ix_] = s2xlib('ii',['F',int2str(round(v_('I')))],ix_);
            pb.xnames{iv} = ['F',int2str(round(v_('I')))];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('nuobservs')
            [ig,ig_] = s2xlib('ii',['RX',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['RX',int2str(I)];
            [ig,ig_] = s2xlib('ii',['RY',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['RY',int2str(I)];
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('RX1')) = -332.65;
        pbm.gconst(ig_('RY1')) = 262.09;
        pbm.gconst(ig_('RX2')) = -199.76;
        pbm.gconst(ig_('RY2')) = 166.7;
        pbm.gconst(ig_('RX3')) = -253.06;
        pbm.gconst(ig_('RY3')) = 202.27;
        pbm.gconst(ig_('RX4')) = 58.13;
        pbm.gconst(ig_('RY4')) = 271.89;
        pbm.gconst(ig_('RX5')) = 238.22;
        pbm.gconst(ig_('RY5')) = 237.37;
        pbm.gconst(ig_('RX6')) = 317.55;
        pbm.gconst(ig_('RY6')) = 221.15;
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = -.6120001572;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = -.6120001572;
        end
        if(isKey(ix_,'Y1'))
            pb.x0(ix_('Y1'),1) = .57175904776;
        else
            pb.y0(find(pbm.congrps==ig_('Y1')),1) = .57175904776;
        end
        if(isKey(ix_,'Z1'))
            pb.x0(ix_('Z1'),1) = -1.847081276;
        else
            pb.y0(find(pbm.congrps==ig_('Z1')),1) = -1.847081276;
        end
        if(isKey(ix_,'RX1'))
            pb.x0(ix_('RX1'),1) = .01574151594;
        else
            pb.y0(find(pbm.congrps==ig_('RX1')),1) = .01574151594;
        end
        if(isKey(ix_,'RY1'))
            pb.x0(ix_('RY1'),1) = -.0127909362;
        else
            pb.y0(find(pbm.congrps==ig_('RY1')),1) = -.0127909362;
        end
        if(isKey(ix_,'RZ1'))
            pb.x0(ix_('RZ1'),1) = -.0044008498;
        else
            pb.y0(find(pbm.congrps==ig_('RZ1')),1) = -.0044008498;
        end
        if(isKey(ix_,'TX1'))
            pb.x0(ix_('TX1'),1) = -.0340938396;
        else
            pb.y0(find(pbm.congrps==ig_('TX1')),1) = -.0340938396;
        end
        if(isKey(ix_,'TY1'))
            pb.x0(ix_('TY1'),1) = -.107513871;
        else
            pb.y0(find(pbm.congrps==ig_('TY1')),1) = -.107513871;
        end
        if(isKey(ix_,'TZ1'))
            pb.x0(ix_('TZ1'),1) = 1.1202240291;
        else
            pb.y0(find(pbm.congrps==ig_('TZ1')),1) = 1.1202240291;
        end
        if(isKey(ix_,'KA1'))
            pb.x0(ix_('KA1'),1) = -3.177064E-7;
        else
            pb.y0(find(pbm.congrps==ig_('KA1')),1) = -3.177064E-7;
        end
        if(isKey(ix_,'KB1'))
            pb.x0(ix_('KB1'),1) = 5.882049E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB1')),1) = 5.882049E-13;
        end
        if(isKey(ix_,'F1'))
            pb.x0(ix_('F1'),1) = 399.75152639;
        else
            pb.y0(find(pbm.congrps==ig_('F1')),1) = 399.75152639;
        end
        if(isKey(ix_,'RX2'))
            pb.x0(ix_('RX2'),1) = .01597732412;
        else
            pb.y0(find(pbm.congrps==ig_('RX2')),1) = .01597732412;
        end
        if(isKey(ix_,'RY2'))
            pb.x0(ix_('RY2'),1) = -.0252244646;
        else
            pb.y0(find(pbm.congrps==ig_('RY2')),1) = -.0252244646;
        end
        if(isKey(ix_,'RZ2'))
            pb.x0(ix_('RZ2'),1) = -.0094001416;
        else
            pb.y0(find(pbm.congrps==ig_('RZ2')),1) = -.0094001416;
        end
        if(isKey(ix_,'TX2'))
            pb.x0(ix_('TX2'),1) = -.0085667661;
        else
            pb.y0(find(pbm.congrps==ig_('TX2')),1) = -.0085667661;
        end
        if(isKey(ix_,'TY2'))
            pb.x0(ix_('TY2'),1) = -.1218804907;
        else
            pb.y0(find(pbm.congrps==ig_('TY2')),1) = -.1218804907;
        end
        if(isKey(ix_,'TZ2'))
            pb.x0(ix_('TZ2'),1) = .7190133075;
        else
            pb.y0(find(pbm.congrps==ig_('TZ2')),1) = .7190133075;
        end
        if(isKey(ix_,'KA2'))
            pb.x0(ix_('KA2'),1) = -3.780477E-7;
        else
            pb.y0(find(pbm.congrps==ig_('KA2')),1) = -3.780477E-7;
        end
        if(isKey(ix_,'KB2'))
            pb.x0(ix_('KB2'),1) = 9.307431E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB2')),1) = 9.307431E-13;
        end
        if(isKey(ix_,'F2'))
            pb.x0(ix_('F2'),1) = 402.01753386;
        else
            pb.y0(find(pbm.congrps==ig_('F2')),1) = 402.01753386;
        end
        if(isKey(ix_,'RX4'))
            pb.x0(ix_('RX4'),1) = .01484625118;
        else
            pb.y0(find(pbm.congrps==ig_('RX4')),1) = .01484625118;
        end
        if(isKey(ix_,'RY4'))
            pb.x0(ix_('RY4'),1) = -.0210628994;
        else
            pb.y0(find(pbm.congrps==ig_('RY4')),1) = -.0210628994;
        end
        if(isKey(ix_,'RZ4'))
            pb.x0(ix_('RZ4'),1) = -.001166948;
        else
            pb.y0(find(pbm.congrps==ig_('RZ4')),1) = -.001166948;
        end
        if(isKey(ix_,'TX4'))
            pb.x0(ix_('TX4'),1) = -.0249509707;
        else
            pb.y0(find(pbm.congrps==ig_('TX4')),1) = -.0249509707;
        end
        if(isKey(ix_,'TY4'))
            pb.x0(ix_('TY4'),1) = -.1139847055;
        else
            pb.y0(find(pbm.congrps==ig_('TY4')),1) = -.1139847055;
        end
        if(isKey(ix_,'TZ4'))
            pb.x0(ix_('TZ4'),1) = .92166020737;
        else
            pb.y0(find(pbm.congrps==ig_('TZ4')),1) = .92166020737;
        end
        if(isKey(ix_,'KA4'))
            pb.x0(ix_('KA4'),1) = -3.295265E-7;
        else
            pb.y0(find(pbm.congrps==ig_('KA4')),1) = -3.295265E-7;
        end
        if(isKey(ix_,'KB4'))
            pb.x0(ix_('KB4'),1) = 6.732885E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB4')),1) = 6.732885E-13;
        end
        if(isKey(ix_,'F4'))
            pb.x0(ix_('F4'),1) = 400.40175368;
        else
            pb.y0(find(pbm.congrps==ig_('F4')),1) = 400.40175368;
        end
        if(isKey(ix_,'RX27'))
            pb.x0(ix_('RX27'),1) = .01991666998;
        else
            pb.y0(find(pbm.congrps==ig_('RX27')),1) = .01991666998;
        end
        if(isKey(ix_,'RY27'))
            pb.x0(ix_('RY27'),1) = -1.22433082;
        else
            pb.y0(find(pbm.congrps==ig_('RY27')),1) = -1.22433082;
        end
        if(isKey(ix_,'RZ27'))
            pb.x0(ix_('RZ27'),1) = .0119988756;
        else
            pb.y0(find(pbm.congrps==ig_('RZ27')),1) = .0119988756;
        end
        if(isKey(ix_,'TX27'))
            pb.x0(ix_('TX27'),1) = -1.411897512;
        else
            pb.y0(find(pbm.congrps==ig_('TX27')),1) = -1.411897512;
        end
        if(isKey(ix_,'TY27'))
            pb.x0(ix_('TY27'),1) = -.1148065151;
        else
            pb.y0(find(pbm.congrps==ig_('TY27')),1) = -.1148065151;
        end
        if(isKey(ix_,'TZ27'))
            pb.x0(ix_('TZ27'),1) = .44915582738;
        else
            pb.y0(find(pbm.congrps==ig_('TZ27')),1) = .44915582738;
        end
        if(isKey(ix_,'KA27'))
            pb.x0(ix_('KA27'),1) = 5.95875E-8;
        else
            pb.y0(find(pbm.congrps==ig_('KA27')),1) = 5.95875E-8;
        end
        if(isKey(ix_,'KB27'))
            pb.x0(ix_('KB27'),1) = -2.48391E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB27')),1) = -2.48391E-13;
        end
        if(isKey(ix_,'F27'))
            pb.x0(ix_('F27'),1) = 407.03024568;
        else
            pb.y0(find(pbm.congrps==ig_('F27')),1) = 407.03024568;
        end
        if(isKey(ix_,'RX30'))
            pb.x0(ix_('RX30'),1) = .02082242153;
        else
            pb.y0(find(pbm.congrps==ig_('RX30')),1) = .02082242153;
        end
        if(isKey(ix_,'RY30'))
            pb.x0(ix_('RY30'),1) = -1.238434791;
        else
            pb.y0(find(pbm.congrps==ig_('RY30')),1) = -1.238434791;
        end
        if(isKey(ix_,'RZ30'))
            pb.x0(ix_('RZ30'),1) = .01389314763;
        else
            pb.y0(find(pbm.congrps==ig_('RZ30')),1) = .01389314763;
        end
        if(isKey(ix_,'TX30'))
            pb.x0(ix_('TX30'),1) = -1.049686225;
        else
            pb.y0(find(pbm.congrps==ig_('TX30')),1) = -1.049686225;
        end
        if(isKey(ix_,'TY30'))
            pb.x0(ix_('TY30'),1) = -.1299513286;
        else
            pb.y0(find(pbm.congrps==ig_('TY30')),1) = -.1299513286;
        end
        if(isKey(ix_,'TZ30'))
            pb.x0(ix_('TZ30'),1) = .33798380231;
        else
            pb.y0(find(pbm.congrps==ig_('TZ30')),1) = .33798380231;
        end
        if(isKey(ix_,'KA30'))
            pb.x0(ix_('KA30'),1) = 4.5673127E-8;
        else
            pb.y0(find(pbm.congrps==ig_('KA30')),1) = 4.5673127E-8;
        end
        if(isKey(ix_,'KB30'))
            pb.x0(ix_('KB30'),1) = -1.79243E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB30')),1) = -1.79243E-13;
        end
        if(isKey(ix_,'F30'))
            pb.x0(ix_('F30'),1) = 405.91764962;
        else
            pb.y0(find(pbm.congrps==ig_('F30')),1) = 405.91764962;
        end
        if(isKey(ix_,'RX37'))
            pb.x0(ix_('RX37'),1) = .01658816461;
        else
            pb.y0(find(pbm.congrps==ig_('RX37')),1) = .01658816461;
        end
        if(isKey(ix_,'RY37'))
            pb.x0(ix_('RY37'),1) = -1.247226838;
        else
            pb.y0(find(pbm.congrps==ig_('RY37')),1) = -1.247226838;
        end
        if(isKey(ix_,'RZ37'))
            pb.x0(ix_('RZ37'),1) = .01846788123;
        else
            pb.y0(find(pbm.congrps==ig_('RZ37')),1) = .01846788123;
        end
        if(isKey(ix_,'TX37'))
            pb.x0(ix_('TX37'),1) = -.8617315756;
        else
            pb.y0(find(pbm.congrps==ig_('TX37')),1) = -.8617315756;
        end
        if(isKey(ix_,'TY37'))
            pb.x0(ix_('TY37'),1) = -.1321089362;
        else
            pb.y0(find(pbm.congrps==ig_('TY37')),1) = -.1321089362;
        end
        if(isKey(ix_,'TZ37'))
            pb.x0(ix_('TZ37'),1) = .28256800868;
        else
            pb.y0(find(pbm.congrps==ig_('TZ37')),1) = .28256800868;
        end
        if(isKey(ix_,'KA37'))
            pb.x0(ix_('KA37'),1) = 4.7465711E-8;
        else
            pb.y0(find(pbm.congrps==ig_('KA37')),1) = 4.7465711E-8;
        end
        if(isKey(ix_,'KB37'))
            pb.x0(ix_('KB37'),1) = -1.50881E-13;
        else
            pb.y0(find(pbm.congrps==ig_('KB37')),1) = -1.50881E-13;
        end
        if(isKey(ix_,'F37'))
            pb.x0(ix_('F37'),1) = 404.73590637;
        else
            pb.y0(find(pbm.congrps==ig_('F37')),1) = 404.73590637;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'eE',iet_);
        elftv{it}{1} = 'RX';
        elftv{it}{2} = 'RY';
        elftv{it}{3} = 'RZ';
        elftv{it}{4} = 'X';
        elftv{it}{5} = 'Y';
        elftv{it}{6} = 'Z';
        elftv{it}{7} = 'TX';
        elftv{it}{8} = 'TY';
        elftv{it}{9} = 'TZ';
        elftv{it}{10} = 'KA';
        elftv{it}{11} = 'KB';
        elftv{it}{12} = 'F';
        elftp{it}{1} = 'YRES';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'EX1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EX2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EX3';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY3';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EX4';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY4';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F27';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EX5';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY5';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F30';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'EX6';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        ename = 'EY6';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eE';
        ielftype(ie) = iet_('eE');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Y1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'Z1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Z',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RX37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RY37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'RZ37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('RZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TX37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TX',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TY37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TY',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'TZ37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('TZ',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KA37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KA',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'KB37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('KB',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'F37';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('F',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('YRES',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('nuobservs')
            ig = ig_(['RX',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EX',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['RY',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EY',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-57-12';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
