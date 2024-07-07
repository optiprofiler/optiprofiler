function varargout = HS103(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Source: problem 103 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: N. Gould, December 1989.
% 
%    classification = 'OOR2-AN-7-5'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS103';

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
        v_('M') = 5;
        v_('N') = 7;
        v_('A101') = 0.5;
        v_('A102') = 0.125;
        v_('A103') = 0.5;
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
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['CONSTR',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['CONSTR',int2str(I)];
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
        pbm.gconst(ig_('CONSTR1')) = 1.0;
        pbm.gconst(ig_('CONSTR2')) = 1.0;
        pbm.gconst(ig_('CONSTR3')) = 1.0;
        pbm.gconst(ig_('CONSTR4')) = 1.0;
        pbm.gconst(ig_('CONSTR5')) = 3000.0;
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.1*ones(pb.n,1);
        pb.xupper = 10.0*ones(pb.n,1);
        pb.xlower(ix_('X7'),1) = 0.01;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 6.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en3PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        [it,iet_] = s2mpjlib( 'ii', 'en4PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        [it,iet_] = s2mpjlib( 'ii', 'en5PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        elftp{it}{5} = 'P5';
        [it,iet_] = s2mpjlib( 'ii', 'en6PR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P6';
        elftp{it}{3} = 'P2';
        elftp{it}{4} = 'P3';
        elftp{it}{5} = 'P4';
        elftp{it}{6} = 'P5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1C1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en4PR';
        ielftype(ie) = iet_('en4PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E2C1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        ename = 'E3C1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.5;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.66666666;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.25;
        ename = 'E1C2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.5;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E2C2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en4PR';
        ielftype(ie) = iet_('en4PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E3C2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.3333333333;
        ename = 'E1C3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.5;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.3333333333;
        ename = 'E2C3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.5;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.5;
        ename = 'E3C3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en4PR';
        ielftype(ie) = iet_('en4PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E4C3';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E1C4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.3333333333;
        ename = 'E2C4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en6PR';
        ielftype(ie) = iet_('en6PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.3333333333;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.666666666;
        [~,posep] = ismember('P6',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.25;
        ename = 'E3C4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -3.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.75;
        ename = 'E4C4';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en3PR';
        ielftype(ie) = iet_('en3PR');
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        ename = 'E1C5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -3.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('A101');
        ename = 'E2C5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en6PR';
        ielftype(ie) = iet_('en6PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P6',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -0.5;
        ename = 'E3C5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en5PR';
        ielftype(ie) = iet_('en5PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        ename = 'E4C5';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'en6PR';
        ielftype(ie) = iet_('en6PR');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,10.0,6.0);
        posev = find(strcmp('V6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.0;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.5;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        [~,posep] = ismember('P6',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 10.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C5');
        pbm.grelw{ig}(posel) = 15.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 20.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4C5');
        pbm.grelw{ig}(posel) = 25.0;
        ig = ig_('CONSTR1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.5;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C1');
        pbm.grelw{ig}(posel) = 0.7;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.2;
        ig = ig_('CONSTR2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.3;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C2');
        pbm.grelw{ig}(posel) = 0.8;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 3.1;
        ig = ig_('CONSTR3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C3');
        pbm.grelw{ig}(posel) = 0.1;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4C3');
        pbm.grelw{ig}(posel) = 0.65;
        ig = ig_('CONSTR4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.2;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C4');
        pbm.grelw{ig}(posel) = 0.3;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 0.4;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4C4');
        pbm.grelw{ig}(posel) = 0.5;
        ig = ig_('CONSTR5');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1C5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 10.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2C5');
        pbm.grelw{ig}(posel) = 15.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3C5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 20.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4C5');
        pbm.grelw{ig}(posel) = 25.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               1809.76476
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = grange(legrps);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-7-5';
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

    case 'en3PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVALUE = (EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^pbm.elpar{iel_}(2))*...
             (EV_(3)^pbm.elpar{iel_}(3));
        varargout{1} = FVALUE;
        if(nargout>1)
            g_(1,1) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1));
            g_(2,1) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2));
            g_(3,1) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) =...
                      FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*((pbm.elpar{iel_}(1)-1.0)/EV_(1));
                H_(2,2) =...
                      FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*((pbm.elpar{iel_}(2)-1.0)/EV_(2));
                H_(3,3) =...
                      FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*((pbm.elpar{iel_}(3)-1.0)/EV_(3));
                H_(1,2) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(2)/EV_(2));
                H_(2,1) = H_(1,2);
                H_(1,3) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,1) = H_(1,3);
                H_(2,3) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    case 'en4PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVALUE = (EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^pbm.elpar{iel_}(2))*...
             (EV_(3)^pbm.elpar{iel_}(3))*(EV_(4)^pbm.elpar{iel_}(4));
        varargout{1} = FVALUE;
        if(nargout>1)
            g_(1,1) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1));
            g_(2,1) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2));
            g_(3,1) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3));
            g_(4,1) = FVALUE*(pbm.elpar{iel_}(4)/EV_(4));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,1) =...
                      FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*((pbm.elpar{iel_}(1)-1.0)/EV_(1));
                H_(2,2) =...
                      FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*((pbm.elpar{iel_}(2)-1.0)/EV_(2));
                H_(3,3) =...
                      FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*((pbm.elpar{iel_}(3)-1.0)/EV_(3));
                H_(4,4) =...
                      FVALUE*(pbm.elpar{iel_}(4)/EV_(4))*((pbm.elpar{iel_}(4)-1.0)/EV_(4));
                H_(1,2) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(2)/EV_(2));
                H_(2,1) = H_(1,2);
                H_(1,3) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,1) = H_(1,3);
                H_(1,4) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,1) = H_(1,4);
                H_(2,3) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,2) = H_(2,3);
                H_(2,4) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,2) = H_(2,4);
                H_(3,4) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,3) = H_(3,4);
                varargout{3} = H_;
            end
        end

    case 'en5PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVALUE = (EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^pbm.elpar{iel_}(2))*...
             (EV_(3)^pbm.elpar{iel_}(3))*(EV_(4)^pbm.elpar{iel_}(4))*(EV_(5)^pbm.elpar{iel_}(5));
        varargout{1} = FVALUE;
        if(nargout>1)
            g_(1,1) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1));
            g_(2,1) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2));
            g_(3,1) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3));
            g_(4,1) = FVALUE*(pbm.elpar{iel_}(4)/EV_(4));
            g_(5,1) = FVALUE*(pbm.elpar{iel_}(5)/EV_(5));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) =...
                      FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*((pbm.elpar{iel_}(1)-1.0)/EV_(1));
                H_(2,2) =...
                      FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*((pbm.elpar{iel_}(2)-1.0)/EV_(2));
                H_(3,3) =...
                      FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*((pbm.elpar{iel_}(3)-1.0)/EV_(3));
                H_(4,4) =...
                      FVALUE*(pbm.elpar{iel_}(4)/EV_(4))*((pbm.elpar{iel_}(4)-1.0)/EV_(4));
                H_(5,5) =...
                      FVALUE*(pbm.elpar{iel_}(5)/EV_(5))*((pbm.elpar{iel_}(5)-1.0)/EV_(5));
                H_(1,2) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(2)/EV_(2));
                H_(2,1) = H_(1,2);
                H_(1,3) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,1) = H_(1,3);
                H_(1,4) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,1) = H_(1,4);
                H_(1,5) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(5)/EV_(5));
                H_(5,1) = H_(1,5);
                H_(2,3) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(3)/EV_(3));
                H_(3,2) = H_(2,3);
                H_(2,4) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,2) = H_(2,4);
                H_(2,5) = FVALUE*(pbm.elpar{iel_}(2)/EV_(2))*(pbm.elpar{iel_}(5)/EV_(5));
                H_(5,2) = H_(2,5);
                H_(3,4) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*(pbm.elpar{iel_}(4)/EV_(4));
                H_(4,3) = H_(3,4);
                H_(3,5) = FVALUE*(pbm.elpar{iel_}(3)/EV_(3))*(pbm.elpar{iel_}(5)/EV_(5));
                H_(5,3) = H_(3,5);
                H_(4,5) = FVALUE*(pbm.elpar{iel_}(4)/EV_(4))*(pbm.elpar{iel_}(5)/EV_(5));
                H_(5,4) = H_(4,5);
                varargout{3} = H_;
            end
        end

    case 'en6PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVALUE = (EV_(1)^pbm.elpar{iel_}(1))*(EV_(2)^pbm.elpar{iel_}(3))*...
             (EV_(3)^pbm.elpar{iel_}(4))*(EV_(4)^pbm.elpar{iel_}(5))*(EV_(5)^pbm.elpar{iel_}(6))*(EV_(6)^pbm.elpar{iel_}(2));
        varargout{1} = FVALUE;
        if(nargout>1)
            g_(1,1) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1));
            g_(2,1) = FVALUE*(pbm.elpar{iel_}(3)/EV_(2));
            g_(3,1) = FVALUE*(pbm.elpar{iel_}(4)/EV_(3));
            g_(4,1) = FVALUE*(pbm.elpar{iel_}(5)/EV_(4));
            g_(5,1) = FVALUE*(pbm.elpar{iel_}(6)/EV_(5));
            g_(6,1) = FVALUE*(pbm.elpar{iel_}(2)/EV_(6));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(6,6);
                H_(1,1) =...
                      FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*((pbm.elpar{iel_}(1)-1.0)/EV_(1));
                H_(2,2) =...
                      FVALUE*(pbm.elpar{iel_}(3)/EV_(2))*((pbm.elpar{iel_}(3)-1.0)/EV_(2));
                H_(3,3) =...
                      FVALUE*(pbm.elpar{iel_}(4)/EV_(3))*((pbm.elpar{iel_}(4)-1.0)/EV_(3));
                H_(4,4) =...
                      FVALUE*(pbm.elpar{iel_}(5)/EV_(4))*((pbm.elpar{iel_}(5)-1.0)/EV_(4));
                H_(5,5) =...
                      FVALUE*(pbm.elpar{iel_}(6)/EV_(5))*((pbm.elpar{iel_}(6)-1.0)/EV_(5));
                H_(6,6) =...
                      FVALUE*(pbm.elpar{iel_}(2)/EV_(6))*((pbm.elpar{iel_}(2)-1.0)/EV_(6));
                H_(1,2) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(3)/EV_(2));
                H_(2,1) = H_(1,2);
                H_(1,3) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(4)/EV_(3));
                H_(3,1) = H_(1,3);
                H_(1,4) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(5)/EV_(4));
                H_(4,1) = H_(1,4);
                H_(1,5) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(6)/EV_(5));
                H_(5,1) = H_(1,5);
                H_(1,6) = FVALUE*(pbm.elpar{iel_}(1)/EV_(1))*(pbm.elpar{iel_}(2)/EV_(6));
                H_(6,1) = H_(1,6);
                H_(2,3) = FVALUE*(pbm.elpar{iel_}(3)/EV_(2))*(pbm.elpar{iel_}(4)/EV_(3));
                H_(3,2) = H_(2,3);
                H_(2,4) = FVALUE*(pbm.elpar{iel_}(3)/EV_(2))*(pbm.elpar{iel_}(5)/EV_(4));
                H_(4,2) = H_(2,4);
                H_(2,5) = FVALUE*(pbm.elpar{iel_}(3)/EV_(2))*(pbm.elpar{iel_}(6)/EV_(5));
                H_(5,2) = H_(2,5);
                H_(2,6) = FVALUE*(pbm.elpar{iel_}(3)/EV_(2))*(pbm.elpar{iel_}(2)/EV_(6));
                H_(6,2) = H_(2,6);
                H_(3,4) = FVALUE*(pbm.elpar{iel_}(4)/EV_(3))*(pbm.elpar{iel_}(5)/EV_(4));
                H_(4,3) = H_(3,4);
                H_(3,5) = FVALUE*(pbm.elpar{iel_}(4)/EV_(3))*(pbm.elpar{iel_}(6)/EV_(5));
                H_(5,3) = H_(3,5);
                H_(3,6) = FVALUE*(pbm.elpar{iel_}(4)/EV_(3))*(pbm.elpar{iel_}(2)/EV_(6));
                H_(6,3) = H_(3,6);
                H_(4,5) = FVALUE*(pbm.elpar{iel_}(5)/EV_(4))*(pbm.elpar{iel_}(6)/EV_(5));
                H_(5,4) = H_(4,5);
                H_(4,6) = FVALUE*(pbm.elpar{iel_}(5)/EV_(4))*(pbm.elpar{iel_}(2)/EV_(6));
                H_(6,4) = H_(4,6);
                H_(5,6) = FVALUE*(pbm.elpar{iel_}(6)/EV_(5))*(pbm.elpar{iel_}(2)/EV_(6));
                H_(6,5) = H_(5,6);
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

