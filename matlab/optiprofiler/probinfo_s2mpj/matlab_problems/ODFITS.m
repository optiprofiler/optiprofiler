function varargout = ODFITS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A simple Origin/Destination matrix fit using a minimum entropy
%    approach.  The objective is a combination of different aims, namely
%    to be close to an a priori matrix for some entries, to be consistent
%    with some traffic counts (for some entries) and to be small (for entries
%    where nothing else is known).
% 
%    The objective function is of the form
%         SUM   m T [ ln( T / a ) - 1 ] + E   SUM  T [ ln ( T  ) - 1 ]
%        i in I  i i       i   i            i in J  i        i
%                +  g   SUM   q  F [ ln( F / c ) - 1 ]
%                     i in K   i  i       i   i
%    with the constraints that all Ti and Fi be positive and that
%                         F  =  SUM p   T
%                          i     j   ij  j
%    where the pij represent path weights from an a priori assignment.
% 
%    Source: a modification of an example in
%    L.G. Willumsen,
%    "Origin-Destination Matrix: static estimation"
%    in "Concise Encyclopedia of Traffic and Transportation Systems"
%    (M. Papageorgiou, ed.), Pergamon Press, 1991.
% 
%    M. Bierlaire, private communication, 1991.
% 
%    SIF input: Ph Toint, Dec 1991.
% 
%    classification = 'C-COLR2-MN-10-6'
% 
%    Number of available traffic counts
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ODFITS';

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
        v_('ARCS') = 6;
        v_('TC1') = 100.0;
        v_('TC2') = 500.0;
        v_('TC3') = 400.0;
        v_('TC4') = 1100.0;
        v_('TC5') = 600.0;
        v_('TC6') = 700.0;
        v_('QLT1') = 1.0;
        v_('QLT2') = 1.0;
        v_('QLT3') = 1.0;
        v_('QLT4') = 1.0;
        v_('QLT5') = 1.0;
        v_('QLT6') = 1.0;
        v_('P131') = 1.0;
        v_('P132') = 0.0;
        v_('P133') = 0.0;
        v_('P134') = 0.0;
        v_('P135') = 0.0;
        v_('P136') = 0.0;
        v_('P141') = 0.0;
        v_('P142') = 1.0;
        v_('P143') = 0.0;
        v_('P144') = 1.0;
        v_('P145') = 0.0;
        v_('P146') = 0.0;
        v_('P231') = 0.0;
        v_('P232') = 0.0;
        v_('P233') = 1.0;
        v_('P234') = 1.0;
        v_('P235') = 1.0;
        v_('P236') = 0.0;
        v_('P241') = 0.0;
        v_('P242') = 0.0;
        v_('P243') = 0.0;
        v_('P244') = 1.0;
        v_('P245') = 1.0;
        v_('P246') = 1.0;
        v_('APV13') = 90.0;
        v_('APV14') = 450.0;
        v_('APV23') = 360.0;
        v_('MU13') = 0.5;
        v_('MU14') = 0.5;
        v_('MU23') = 0.5;
        v_('1/MU13') = 1.0/v_('MU13');
        v_('1/MU14') = 1.0/v_('MU14');
        v_('1/MU23') = 1.0/v_('MU23');
        v_('GAMMA') = 1.5;
        v_('ENTROP') = 0.2;
        v_('1/ENTR') = 1.0/v_('ENTROP');
        v_('1') = 1;
        for I=v_('1'):v_('ARCS')
            v_(['1/QLT',int2str(I)]) = 1.0/v_(['QLT',int2str(I)]);
            v_(['G/QLT',int2str(I)]) = v_(['1/QLT',int2str(I)])*v_('GAMMA');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','T13',ix_);
        pb.xnames{iv} = 'T13';
        [iv,ix_] = s2mpjlib('ii','T14',ix_);
        pb.xnames{iv} = 'T14';
        [iv,ix_] = s2mpjlib('ii','T23',ix_);
        pb.xnames{iv} = 'T23';
        [iv,ix_] = s2mpjlib('ii','T24',ix_);
        pb.xnames{iv} = 'T24';
        for I=v_('1'):v_('ARCS')
            [iv,ix_] = s2mpjlib('ii',['F',int2str(I)],ix_);
            pb.xnames{iv} = ['F',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','AP13',ig_);
        gtype{ig} = '<>';
        iv = ix_('T13');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('1/MU13');
        [ig,ig_] = s2mpjlib('ii','AP14',ig_);
        gtype{ig} = '<>';
        iv = ix_('T14');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('1/MU14');
        [ig,ig_] = s2mpjlib('ii','AP23',ig_);
        gtype{ig} = '<>';
        iv = ix_('T23');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('1/MU23');
        [ig,ig_] = s2mpjlib('ii','AP24',ig_);
        gtype{ig} = '<>';
        iv = ix_('T24');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        pbm.gscale(ig,1) = v_('1/ENTR');
        for I=v_('1'):v_('ARCS')
            [ig,ig_] = s2mpjlib('ii',['CP',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['F',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            pbm.gscale(ig,1) = v_(['G/QLT',int2str(I)]);
        end
        for I=v_('1'):v_('ARCS')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(I)];
            iv = ix_(['F',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_('T13');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['P13',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['P13',int2str(I)]);
            end
            iv = ix_('T14');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['P14',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['P14',int2str(I)]);
            end
            iv = ix_('T23');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['P23',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['P23',int2str(I)]);
            end
            iv = ix_('T24');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_(['P24',int2str(I)])+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_(['P24',int2str(I)]);
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.1*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('T13'),1) = v_('APV13');
        pb.x0(ix_('T14'),1) = v_('APV14');
        pb.x0(ix_('T23'),1) = v_('APV23');
        pb.x0(ix_('T24'),1) = 1.0;
        for I=v_('1'):v_('ARCS')
            pb.x0(ix_(['F',int2str(I)]),1) = v_(['TC',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eXLOGX',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'DEN';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'TFIT13';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXLOGX';
        ielftype(ie) = iet_('eXLOGX');
        vname = 'T13';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('DEN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('APV13');
        ename = 'TFIT23';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXLOGX';
        ielftype(ie) = iet_('eXLOGX');
        vname = 'T23';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('DEN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('APV23');
        ename = 'TFIT14';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXLOGX';
        ielftype(ie) = iet_('eXLOGX');
        vname = 'T14';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('DEN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('APV14');
        ename = 'TFIT24';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eXLOGX';
        ielftype(ie) = iet_('eXLOGX');
        vname = 'T24';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('DEN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1.0;
        for I=v_('1'):v_('ARCS')
            ename = ['CFIT',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eXLOGX';
            ielftype(ie) = iet_('eXLOGX');
            vname = ['F',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.1,[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('DEN',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['TC',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('AP13');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TFIT13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('AP14');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TFIT14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('AP23');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TFIT23');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('AP24');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('TFIT24');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('1'):v_('ARCS')
            ig = ig_(['CP',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['CFIT',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO ODFITS             -2380.026775
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-MN-10-6';
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

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eXLOGX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LOGX = log(EV_(1)/pbm.elpar{iel_}(1));
        varargout{1} = EV_(1)*LOGX;
        if(nargout>1)
            g_(1,1) = 1.0+LOGX;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0/EV_(1);
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

