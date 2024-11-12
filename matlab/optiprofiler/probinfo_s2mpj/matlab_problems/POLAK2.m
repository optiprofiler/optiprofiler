function varargout = POLAK2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : POLAK2
%    *********
% 
%    A nonlinear minmax problem in ten variables.
% 
%    Source: 
%    E. Polak, D.H. Mayne and J.E. Higgins,
%    "Superlinearly convergent algorithm for min-max problems"
%    JOTA 69, pp. 407-439, 1991.
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'C-CLOR2-AN-11-2'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'POLAK2';

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
        v_('1') = 1;
        v_('10') = 10;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('10')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii','F1',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F1';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii','F2',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'F2';
        iv = ix_('U');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
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
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.1*ones(pb.n,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 100.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 100.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEL',iet_);
        elftv{it}{1} = 'XX1';
        elftv{it}{2} = 'XX2';
        elftv{it}{3} = 'XX3';
        elftv{it}{4} = 'XX4';
        elftv{it}{5} = 'XX5';
        elftv{it}{6} = 'XX6';
        elftv{it}{7} = 'XX7';
        elftv{it}{8} = 'XX8';
        elftv{it}{9} = 'XX9';
        elftv{it}{10} = 'XX10';
        elftp{it}{1} = 'P';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL';
        ielftype(ie) = iet_('eEL');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 2.0;
        ename = 'E2';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eEL';
        ielftype(ie) = iet_('eEL');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX4',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X5';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX5',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X6';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX6',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X7';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX7',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X8';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX8',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X9';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX9',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X10';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1);
        posev = find(strcmp('XX10',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.0;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('F1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('F2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               54.598146
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-AN-11-2';
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

    case 'eEL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = 1.0e-8*EV_(1)*EV_(1)+(EV_(2)+pbm.elpar{iel_}(1))^2;
        A = A+EV_(3)*EV_(3)+4.0*EV_(4)*EV_(4);
        A = A+EV_(5)*EV_(5)+EV_(6)*EV_(6)+EV_(7)*EV_(7);
        A = A+EV_(8)*EV_(8)+EV_(9)*EV_(9)+EV_(10)*EV_(10);
        EA = exp(A);
        varargout{1} = EA;
        if(nargout>1)
            g_(1,1) = 2.0e-8*EV_(1)*EA;
            g_(2,1) = 2.0*(EV_(2)+pbm.elpar{iel_}(1))*EA;
            g_(3,1) = 2.0*EV_(3)*EA;
            g_(4,1) = 8.0*EV_(4)*EA;
            g_(5,1) = 2.0*EV_(5)*EA;
            g_(6,1) = 2.0*EV_(6)*EA;
            g_(7,1) = 2.0*EV_(7)*EA;
            g_(8,1) = 2.0*EV_(8)*EA;
            g_(9,1) = 2.0*EV_(9)*EA;
            g_(10,1) = 2.0*EV_(10)*EA;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(10,10);
                H_(1,1) = 2.0e-8*EA*(1.0+2.0e-8*EV_(1)^2);
                H_(1,2) = 4.0e-8*EV_(1)*(EV_(2)+pbm.elpar{iel_}(1))*EA;
                H_(2,1) = H_(1,2);
                H_(1,3) = 4.0e-8*EV_(1)*EV_(3)*EA;
                H_(3,1) = H_(1,3);
                H_(1,4) = 1.6e-7*EV_(1)*EV_(4)*EA;
                H_(4,1) = H_(1,4);
                H_(1,5) = 4.0e-8*EV_(1)*EV_(5)*EA;
                H_(5,1) = H_(1,5);
                H_(1,6) = 4.0e-8*EV_(1)*EV_(6)*EA;
                H_(6,1) = H_(1,6);
                H_(1,7) = 4.0e-8*EV_(1)*EV_(7)*EA;
                H_(7,1) = H_(1,7);
                H_(1,8) = 4.0e-8*EV_(1)*EV_(8)*EA;
                H_(8,1) = H_(1,8);
                H_(1,9) = 4.0e-8*EV_(1)*EV_(9)*EA;
                H_(9,1) = H_(1,9);
                H_(1,10) = 4.0e-8*EV_(1)*EV_(10)*EA;
                H_(10,1) = H_(1,10);
                H_(2,2) = 2.0*EA*(1.0+2.0*(EV_(2)+pbm.elpar{iel_}(1))^2);
                H_(2,3) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(3)*EA;
                H_(3,2) = H_(2,3);
                H_(2,4) = 16.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(4)*EA;
                H_(4,2) = H_(2,4);
                H_(2,5) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(5)*EA;
                H_(5,2) = H_(2,5);
                H_(2,6) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(6)*EA;
                H_(6,2) = H_(2,6);
                H_(2,7) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(7)*EA;
                H_(7,2) = H_(2,7);
                H_(2,8) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(8)*EA;
                H_(8,2) = H_(2,8);
                H_(2,9) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(9)*EA;
                H_(9,2) = H_(2,9);
                H_(2,10) = 4.0*(EV_(2)+pbm.elpar{iel_}(1))*EV_(10)*EA;
                H_(10,2) = H_(2,10);
                H_(3,3) = 2.0*EA*(1.0+2.0*EV_(3)*EV_(3));
                H_(3,4) = 16.0*EV_(3)*EV_(4)*EA;
                H_(4,3) = H_(3,4);
                H_(3,5) = 4.0*EV_(3)*EV_(5)*EA;
                H_(5,3) = H_(3,5);
                H_(3,6) = 4.0*EV_(3)*EV_(6)*EA;
                H_(6,3) = H_(3,6);
                H_(3,7) = 4.0*EV_(3)*EV_(7)*EA;
                H_(7,3) = H_(3,7);
                H_(3,8) = 4.0*EV_(3)*EV_(8)*EA;
                H_(8,3) = H_(3,8);
                H_(3,9) = 4.0*EV_(3)*EV_(9)*EA;
                H_(9,3) = H_(3,9);
                H_(3,10) = 4.0*EV_(3)*EV_(10)*EA;
                H_(10,3) = H_(3,10);
                H_(4,4) = 8.0*EA*(1.0+8.0*EV_(4)*EV_(4));
                H_(4,5) = 16.0*EV_(4)*EV_(5)*EA;
                H_(5,4) = H_(4,5);
                H_(4,6) = 16.0*EV_(4)*EV_(6)*EA;
                H_(6,4) = H_(4,6);
                H_(4,7) = 16.0*EV_(4)*EV_(7)*EA;
                H_(7,4) = H_(4,7);
                H_(4,8) = 16.0*EV_(4)*EV_(8)*EA;
                H_(8,4) = H_(4,8);
                H_(4,9) = 16.0*EV_(4)*EV_(9)*EA;
                H_(9,4) = H_(4,9);
                H_(4,10) = 16.0*EV_(4)*EV_(10)*EA;
                H_(10,4) = H_(4,10);
                H_(5,5) = 2.0*EA*(1.0+2.0*EV_(5)*EV_(5));
                H_(5,6) = 4.0*EV_(5)*EV_(6)*EA;
                H_(6,5) = H_(5,6);
                H_(5,7) = 4.0*EV_(5)*EV_(7)*EA;
                H_(7,5) = H_(5,7);
                H_(5,8) = 4.0*EV_(5)*EV_(8)*EA;
                H_(8,5) = H_(5,8);
                H_(5,9) = 4.0*EV_(5)*EV_(9)*EA;
                H_(9,5) = H_(5,9);
                H_(5,10) = 4.0*EV_(5)*EV_(10)*EA;
                H_(10,5) = H_(5,10);
                H_(6,6) = 2.0*EA*(1.0+2.0*EV_(6)*EV_(6));
                H_(6,7) = 4.0*EV_(6)*EV_(7)*EA;
                H_(7,6) = H_(6,7);
                H_(6,8) = 4.0*EV_(6)*EV_(8)*EA;
                H_(8,6) = H_(6,8);
                H_(6,9) = 4.0*EV_(6)*EV_(9)*EA;
                H_(9,6) = H_(6,9);
                H_(6,10) = 4.0*EV_(6)*EV_(10)*EA;
                H_(10,6) = H_(6,10);
                H_(7,7) = 2.0*EA*(1.0+2.0*EV_(7)*EV_(7));
                H_(7,8) = 4.0*EV_(7)*EV_(8)*EA;
                H_(8,7) = H_(7,8);
                H_(7,9) = 4.0*EV_(7)*EV_(9)*EA;
                H_(9,7) = H_(7,9);
                H_(7,10) = 4.0*EV_(7)*EV_(10)*EA;
                H_(10,7) = H_(7,10);
                H_(8,8) = 2.0*EA*(1.0+2.0*EV_(8)*EV_(8));
                H_(8,9) = 4.0*EV_(8)*EV_(9)*EA;
                H_(9,8) = H_(8,9);
                H_(8,10) = 4.0*EV_(8)*EV_(10)*EA;
                H_(10,8) = H_(8,10);
                H_(9,9) = 2.0*EA*(1.0+2.0*EV_(9)*EV_(9));
                H_(9,10) = 4.0*EV_(9)*EV_(10)*EA;
                H_(10,9) = H_(9,10);
                H_(10,10) = 2.0*EA*(1.0+2.0*EV_(10)*EV_(10));
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

