function varargout = HS69(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS69
%    *********
% 
%    This is a cost optimal inspection plan.
% 
%    Source: problem 69 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'OOR2-MN-4-2'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS69';

switch(action)

    case 'setup'

    pb.name      = 'HS69';
    pb.sifpbname = 'HS69';
    pbm.name     = 'HS69';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('N') = 4;
        v_('A') = 0.1;
        v_('B') = 1000.0;
        v_('D') = 1.0;
        v_('NN') = 4.0;
        v_('1') = 1;
        v_('AN') = v_('A')*v_('NN');
        v_('ROOTN') = sqrt(v_('NN'));
        v_('DROOTN') = v_('D')*v_('ROOTN');
        v_('-DROOTN') = -1.0*v_('DROOTN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2xlib('ii','C1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C1';
        iv = ix_('X3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0e+0;
        end
        [ig,ig_] = s2xlib('ii','C2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'C2';
        iv = ix_('X4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0e+0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0e+0;
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.0001;
        pb.xupper(ix_('X1')) = 100.0;
        pb.xupper(ix_('X2')) = 100.0;
        pb.xupper(ix_('X3')) = 2.0;
        pb.xupper(ix_('X4')) = 2.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.y0 = 1.0*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'RECIP',iet_);
        elftv{it}{1} = 'X1';
        [it,iet_] = s2xlib( 'ii', 'NASTYEXP',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X3';
        elftv{it}{3} = 'X4';
        elftp{it}{1} = 'B';
        [it,iet_] = s2xlib( 'ii', 'PHI',iet_);
        elftv{it}{1} = 'X2';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'OE1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'RECIP';
        ielftype(ie) = iet_('RECIP');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'OE2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'NASTYEXP';
        ielftype(ie) = iet_('NASTYEXP');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('B',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('B');
        ename = 'C1E1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PHI';
        ielftype(ie) = iet_('PHI');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0e+0;
        ename = 'C2E1';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PHI';
        ielftype(ie) = iet_('PHI');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('DROOTN');
        ename = 'C2E2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'PHI';
        ielftype(ie) = iet_('PHI');
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,[],[],1.0);
        posev = find(strcmp('X2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('-DROOTN');
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('AN');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('OE2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0e+0;
        ig = ig_('C1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C1E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -2.0e+0;
        ig = ig_('C2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('C2E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -1.0e+0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('C2E2');
        pbm.grelw{ig}(posel) = -1.0e+0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-4-2';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
