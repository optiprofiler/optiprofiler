function varargout = HS67(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS67
%    *********
% 
%    Source: problem 67 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    Original Source: problem 8 in
%    A.R. Colville
%    "A comparative study on nonlinear programming"
%    IBM Scientific Center Report 320-2949, New York, 1968.
% 
%    SIF input: A.R. Conn & Nick Gould, April 1991.
% 
%    classification = 'OOI2-AN-3-14'
% 
%    Set useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS67';

switch(action)

    case 'setup'

    pb.name      = 'HS67';
    pb.sifpbname = 'HS67';
    pbm.name     = 'HS67';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_(['A',int2str(v_('1'))]) = 0.0;
        v_(['A',int2str(v_('2'))]) = 0.0;
        v_(['A',int2str(v_('3'))]) = 85.0;
        v_(['A',int2str(v_('4'))]) = 90.0;
        v_(['A',int2str(v_('5'))]) = 3.0;
        v_(['A',int2str(v_('6'))]) = 0.01;
        v_(['A',int2str(v_('7'))]) = 145.0;
        v_(['A',int2str(v_('8'))]) = 5000.0;
        v_(['A',int2str(v_('9'))]) = 2000.0;
        v_(['A',int2str(v_('10'))]) = 93.0;
        v_(['A',int2str(v_('11'))]) = 95.0;
        v_(['A',int2str(v_('12'))]) = 12.0;
        v_(['A',int2str(v_('13'))]) = 4.0;
        v_(['A',int2str(v_('14'))]) = 162.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('3')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2xlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(v_('1'))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 5.04+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 5.04;
        end
        iv = ix_(['X',int2str(v_('2'))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 0.035+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 0.035;
        end
        iv = ix_(['X',int2str(v_('3'))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 10.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 10.0;
        end
        for I=v_('1'):v_('7')
            [ig,ig_] = s2xlib('ii',['AG',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['AG',int2str(I)];
            [ig,ig_] = s2xlib('ii',['AL',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['AL',int2str(I)];
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
        for I=v_('1'):v_('7')
            v_('I+7') = 7+I;
            pbm.gconst(ig_(['AG',int2str(I)])) = v_(['A',int2str(I)]);
            pbm.gconst(ig_(['AL',int2str(I)])) = v_(['A',int2str(v_('I+7'))]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.00001*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xupper(ix_(['X',int2str(v_('1'))])) = 2000.0;
        pb.xupper(ix_(['X',int2str(v_('2'))])) = 16000.0;
        pb.xupper(ix_(['X',int2str(v_('3'))])) = 120.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(v_('1'))]),1) = 1745.0;
        pb.x0(ix_(['X',int2str(v_('2'))]),1) = 12000.0;
        pb.x0(ix_(['X',int2str(v_('3'))]),1) = 110.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2xlib( 'ii', 'Y2Y5',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y2',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y3',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y4',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y5',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y6',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y7',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        [it,iet_] = s2xlib( 'ii', 'Y8',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        ename = 'E25';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y2Y5';
        ielftype(ie) = iet_('Y2Y5');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E2';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y2';
        ielftype(ie) = iet_('Y2');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E3';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y3';
        ielftype(ie) = iet_('Y3');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E4';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y4';
        ielftype(ie) = iet_('Y4');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E5';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y5';
        ielftype(ie) = iet_('Y5');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E6';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y6';
        ielftype(ie) = iet_('Y6');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E7';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y7';
        ielftype(ie) = iet_('Y7');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        ename = 'E8';
        [ie,ie_] = s2xlib('ii',ename,ie_);
        pbm.elftype{ie} = 'Y8';
        ielftype(ie) = iet_('Y8');
        vname = 'X1';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2xlib('nlx',vname,ix_,pb,1,0.00001,[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = -0.063;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        pbm.grelw{ig}(posel) = 3.36;
        ig = ig_(['AG',int2str(v_('1'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('1'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('2'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('2'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('3'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('3'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('4'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('4'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('5'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('5'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('6'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('6'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AG',int2str(v_('7'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        ig = ig_(['AL',int2str(v_('7'))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTIONS %%%%%%%
        [~,lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOI2-AN-3-14';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%
