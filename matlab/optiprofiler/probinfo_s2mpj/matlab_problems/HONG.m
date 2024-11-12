function varargout = HONG(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    Source: Se June Hong/Chid Apte
% 
%    SIF input: A.R.Conn, Jan 1991.
% 
%    classification = 'C-COLR2-AN-4-1'
% 
%   Problem parameters
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HONG';

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
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','T1',ix_);
        pb.xnames{iv} = 'T1';
        [iv,ix_] = s2mpjlib('ii','T2',ix_);
        pb.xnames{iv} = 'T2';
        [iv,ix_] = s2mpjlib('ii','T3',ix_);
        pb.xnames{iv} = 'T3';
        [iv,ix_] = s2mpjlib('ii','T4',ix_);
        pb.xnames{iv} = 'T4';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','SUM1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'SUM1';
        iv = ix_('T1');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('T2');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('T3');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        iv = ix_('T4');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('SUM1')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.5*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXP',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'P1';
        elftp{it}{2} = 'P2';
        elftp{it}{3} = 'P3';
        elftp{it}{4} = 'P4';
        elftp{it}{5} = 'P5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E1';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
        end
        vname = 'T1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,0.5);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 25.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.92;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.08;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.38;
        ename = 'E2';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
        end
        vname = 'T2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,0.5);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 50.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -2.95;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 3.95;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.11;
        ename = 'E3';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
        end
        vname = 'T3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,0.5);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -4.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.66;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 1657834.0;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = -1.48;
        ename = 'E4';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
        end
        vname = 'T4';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,0.5);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('P1',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.0;
        [~,posep] = ismember('P2',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 20000.0;
        [~,posep] = ismember('P3',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.11;
        [~,posep] = ismember('P4',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.89;
        [~,posep] = ismember('P5',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 0.00035;
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution unknown
        pb.objlower = -4.0;
        pb.objupper = 300.0;
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AN-4-1';
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

    case 'eEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        XTOT = pbm.elpar{iel_}(1)+pbm.elpar{iel_}(2)*EV_(1);
        EP5 = exp(pbm.elpar{iel_}(5)*XTOT);
        varargout{1} = pbm.elpar{iel_}(3)+pbm.elpar{iel_}(4)*EP5;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(4)*pbm.elpar{iel_}(5)*EP5;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = pbm.elpar{iel_}(2)*pbm.elpar{iel_}(2)*pbm.elpar{iel_}(4)*...
                     pbm.elpar{iel_}(5)*pbm.elpar{iel_}(5)*EP5;
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

