function varargout = ORTHREGB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ORTHREGB
%    *********
% 
%    An orthogonal regression problem.
% 
%    The problem is to fit (orthogonally) an ellipse to a set of 6 points
%    in the 3D space. These points are compatible with this constraint.
% 
%    Source:
%    M. Gulliksson,
%    "Algorithms for nonlinear Least-squares with Applications to
%    Orthogonal Regression",
%    UMINF-178.90, University of Umea, Sweden, 1990.
% 
%    SIF input: Ph. Toint, June 1990.
%               correction by Ph. Shott, Jan 1995.
% 
%    classification = 'QQR2-AN-27-6'
% 
%    Parameters for the generation of the data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ORTHREGB';

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
        v_('A') = 9.0;
        v_('B') = 6.0;
        v_('C') = 7.0;
        v_('CX') = 0.5;
        v_('CY') = 0.5;
        v_('CZ') = 0.5;
        v_('1') = 1;
        v_('-A') = -1.0*v_('A');
        v_('-B') = -1.0*v_('B');
        v_('-C') = -1.0*v_('C');
        v_('NPTS') = 1;
        v_('XZ') = v_('CX');
        v_('YZ') = v_('CY');
        v_('ZZ') = v_('CZ');
        v_('NPTS') = 1;
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ')+v_('A');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ')+v_('A');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ');
        v_('NPTS') = 1+v_('NPTS');
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ')+v_('B');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ')+v_('-B');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ');
        v_('NPTS') = 1+v_('NPTS');
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ')+v_('-A');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ')+v_('-A');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ');
        v_('NPTS') = 1+v_('NPTS');
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ')+v_('-B');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ')+v_('B');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ');
        v_('NPTS') = 1+v_('NPTS');
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ')+v_('C');
        v_('NPTS') = 1+v_('NPTS');
        v_(['XD',int2str(round(v_('NPTS')))]) = v_('XZ');
        v_(['YD',int2str(round(v_('NPTS')))]) = v_('YZ');
        v_(['ZD',int2str(round(v_('NPTS')))]) = v_('ZZ')+v_('-C');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','H11',ix_);
        pb.xnames{iv} = 'H11';
        [iv,ix_] = s2mpjlib('ii','H12',ix_);
        pb.xnames{iv} = 'H12';
        [iv,ix_] = s2mpjlib('ii','H13',ix_);
        pb.xnames{iv} = 'H13';
        [iv,ix_] = s2mpjlib('ii','H22',ix_);
        pb.xnames{iv} = 'H22';
        [iv,ix_] = s2mpjlib('ii','H23',ix_);
        pb.xnames{iv} = 'H23';
        [iv,ix_] = s2mpjlib('ii','H33',ix_);
        pb.xnames{iv} = 'H33';
        [iv,ix_] = s2mpjlib('ii','G1',ix_);
        pb.xnames{iv} = 'G1';
        [iv,ix_] = s2mpjlib('ii','G2',ix_);
        pb.xnames{iv} = 'G2';
        [iv,ix_] = s2mpjlib('ii','G3',ix_);
        pb.xnames{iv} = 'G3';
        for I=v_('1'):v_('NPTS')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NPTS')
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['OY',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['OZ',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['Z',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(I)];
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
        for I=v_('1'):v_('NPTS')
            pbm.gconst(ig_(['OX',int2str(I)])) = v_(['XD',int2str(I)]);
            pbm.gconst(ig_(['OY',int2str(I)])) = v_(['YD',int2str(I)]);
            pbm.gconst(ig_(['OZ',int2str(I)])) = v_(['ZD',int2str(I)]);
            pbm.gconst(ig_(['E',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'H11'))
            pb.x0(ix_('H11'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('H11')),1) = 1.0;
        end
        if(isKey(ix_,'H12'))
            pb.x0(ix_('H12'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('H12')),1) = 0.0;
        end
        if(isKey(ix_,'H13'))
            pb.x0(ix_('H13'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('H13')),1) = 0.0;
        end
        if(isKey(ix_,'H22'))
            pb.x0(ix_('H22'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('H22')),1) = 1.0;
        end
        if(isKey(ix_,'H23'))
            pb.x0(ix_('H23'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('H23')),1) = 0.0;
        end
        if(isKey(ix_,'H33'))
            pb.x0(ix_('H33'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('H33')),1) = 1.0;
        end
        if(isKey(ix_,'G1'))
            pb.x0(ix_('G1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('G1')),1) = 0.0;
        end
        if(isKey(ix_,'G2'))
            pb.x0(ix_('G2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('G2')),1) = 0.0;
        end
        if(isKey(ix_,'G3'))
            pb.x0(ix_('G3'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('G3')),1) = 0.0;
        end
        for I=v_('1'):v_('NPTS')
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_(['XD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_(['XD',int2str(I)]);
            end
            if(isKey(ix_,['Y',int2str(I)]))
                pb.x0(ix_(['Y',int2str(I)]),1) = v_(['YD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['Y',int2str(I)])),1) = v_(['YD',int2str(I)]);
            end
            if(isKey(ix_,['Z',int2str(I)]))
                pb.x0(ix_(['Z',int2str(I)]),1) = v_(['ZD',int2str(I)]);
            else
                pb.y0(find(pbm.congrps==ig_(['Z',int2str(I)])),1) = v_(['ZD',int2str(I)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eHXX',iet_);
        elftv{it}{1} = 'H';
        elftv{it}{2} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eHXY',iet_);
        elftv{it}{1} = 'H';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eGX',iet_);
        elftv{it}{1} = 'G';
        elftv{it}{2} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NPTS')
            ename = ['EA',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXX';
            ielftype(ie) = iet_('eHXX');
            vname = 'H11';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXY';
            ielftype(ie) = iet_('eHXY');
            vname = 'H12';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXX';
            ielftype(ie) = iet_('eHXX');
            vname = 'H22';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['ED',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eGX';
            ielftype(ie) = iet_('eGX');
            vname = 'G1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('G',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EE',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eGX';
            ielftype(ie) = iet_('eGX');
            vname = 'G2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('G',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EF',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXY';
            ielftype(ie) = iet_('eHXY');
            vname = 'H13';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EG',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXY';
            ielftype(ie) = iet_('eHXY');
            vname = 'H23';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EH',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eHXX';
            ielftype(ie) = iet_('eHXX');
            vname = 'H33';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EI',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eGX';
            ielftype(ie) = iet_('eGX');
            vname = 'G3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('G',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NPTS')
            ig = ig_(['OX',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OY',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OZ',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['E',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(I)]);
            pbm.grelw{ig}(posel) = 2.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['ED',int2str(I)]);
            pbm.grelw{ig}(posel) = -2.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EE',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -2.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EF',int2str(I)]);
            pbm.grelw{ig}(posel) = 2.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EG',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 2.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EH',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EI',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -2.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QQR2-AN-27-6';
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

    case 'eHXX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(2);
            g_(2,1) = 2.0*EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = EV_(2)+EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)+EV_(1);
                varargout{3} = H_;
            end
        end

    case 'eHXY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2)*EV_(3);
        if(nargout>1)
            g_(1,1) = EV_(2)*EV_(3);
            g_(2,1) = EV_(1)*EV_(3);
            g_(3,1) = EV_(1)*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = EV_(3);
                H_(2,1) = H_(1,2);
                H_(1,3) = EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,3) = EV_(1);
                H_(3,2) = H_(2,3);
                varargout{3} = H_;
            end
        end

    case 'eGX'

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

