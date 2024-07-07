function varargout = MANNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MANNE
%    *********
% 
%    A variable dimension econometric equilibrium problem
%    suggested by A. Manne
% 
%    Source:
%    B. Murtagh and M. Saunders,
%    Mathematical Programming Studies 16, pp. 84-117,
%    (example 5.12).
% 
%    SIF input: N. Gould and Ph. Toint, March 1990.
% 
%    classification = 'OOR2-MN-V-V'
% 
%    Number of periods
%    The number of variables in the problem N = 3*T
% 
%       Alternative values for the SIF file parameters:
% IE T                   100            $-PARAMETER n = 300    original value
% IE T                   365            $-PARAMETER n = 995
% IE T                   1000           $-PARAMETER n = 3000
% IE T                   2000           $-PARAMETER n = 6000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MANNE';

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
        if(nargs<1)
            v_('T') = 4;  %  SIF file default value
        else
            v_('T') = varargin{1};
        end
        v_('GROW') = 0.03;
        v_('BETA') = 0.95;
        v_('XK0') = 3.0;
        v_('XC0') = 0.95;
        v_('XI0') = 0.05;
        v_('B') = 0.25;
        v_('BPROB') = 1.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('T-1') = -1+v_('T');
        v_('T-2') = -2+v_('T');
        v_('LOGXK') = log(v_('XK0'));
        v_('BLOGX') = v_('LOGXK')*v_('B');
        v_('XK0**B') = exp(v_('BLOGX'));
        v_('NUM') = v_('XC0')+v_('XI0');
        v_('A') = v_('NUM')/v_('XK0**B');
        v_('1-B') = 1.0-v_('B');
        v_('1+G') = 1.0+v_('GROW');
        v_('LOG1+G') = log(v_('1+G'));
        v_('SOME') = v_('LOG1+G')*v_('1-B');
        v_('GFAC') = exp(v_('SOME'));
        v_('AT1') = v_('A')*v_('GFAC');
        v_('BT1') = 0.0+v_('BETA');
        for J=v_('2'):v_('T')
            v_('J-1') = -1+J;
            v_(['AT',int2str(J)]) = v_(['AT',int2str(round(v_('J-1')))])*v_('GFAC');
            v_(['BT',int2str(J)]) = v_(['BT',int2str(round(v_('J-1')))])*v_('BETA');
        end
        v_('1-BETA') = 1.0-v_('BETA');
        v_('1/1-BETA') = 1.0/v_('1-BETA');
        v_(['BT',int2str(round(v_('T')))]) = v_(['BT',int2str(round(v_('T')))])*...
             v_('1/1-BETA');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('T')
            [iv,ix_] = s2mpjlib('ii',['C',int2str(I)],ix_);
            pb.xnames{iv} = ['C',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['I',int2str(I)],ix_);
            pb.xnames{iv} = ['I',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['K',int2str(I)],ix_);
            pb.xnames{iv} = ['K',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('T')
            [ig,ig_] = s2mpjlib('ii',['NL',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['NL',int2str(I)];
            iv = ix_(['C',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['I',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for I=v_('1'):v_('T-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['L',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['L',int2str(I)];
            iv = ix_(['K',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['K',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['I',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii',['L',int2str(round(v_('T')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['L',int2str(round(v_('T')))];
        iv = ix_(['K',int2str(round(v_('T')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('GROW')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('GROW');
        end
        [ig,ig_] = s2mpjlib('ii',['L',int2str(round(v_('T')))],ig_);
        gtype{ig}  = '<=';
        cnames{ig} = ['L',int2str(round(v_('T')))];
        iv = ix_(['I',int2str(round(v_('T')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
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
        %%%%%%%%%%%%%%%%%%%%%  RANGES %%%%%%%%%%%%%%%%%%%%%%
        grange(legrps,1) = Inf*ones(pb.nle,1);
        grange(gegrps,1) = Inf*ones(pb.nge,1);
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('K1'),1) = 3.05;
        pb.xupper(ix_('K1'),1) = 3.05;
        for I=v_('2'):v_('T')
            pb.xlower(ix_(['K',int2str(I)]),1) = 3.05;
        end
        v_('1.04**T') = 0.05;
        for I=v_('1'):v_('T')
            v_('1.04**T') = 1.04*v_('1.04**T');
            pb.xlower(ix_(['C',int2str(I)]),1) = 0.95;
            pb.xlower(ix_(['I',int2str(I)]),1) = 0.05;
            pb.xupper(ix_(['I',int2str(I)])) = v_('1.04**T');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'K1'))
            pb.x0(ix_('K1'),1) = 3.05;
        else
            pb.y0(find(pbm.congrps==ig_('K1')),1) = 3.05;
        end
        for I=v_('2'):v_('T')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('I-1/10') = 0.1*v_('RI-1');
            v_('VAL') = 3.0+v_('I-1/10');
            if(isKey(ix_,['K',int2str(I)]))
                pb.x0(ix_(['K',int2str(I)]),1) = v_('VAL');
            else
                pb.y0(find(pbm.congrps==ig_(['K',int2str(I)])),1) = v_('VAL');
            end
        end
        for I=v_('1'):v_('T')
            if(isKey(ix_,['C',int2str(I)]))
                pb.x0(ix_(['C',int2str(I)]),1) = 0.95;
            else
                pb.y0(find(pbm.congrps==ig_(['C',int2str(I)])),1) = 0.95;
            end
            if(isKey(ix_,['I',int2str(I)]))
                pb.x0(ix_(['I',int2str(I)]),1) = 0.05;
            else
                pb.y0(find(pbm.congrps==ig_(['I',int2str(I)])),1) = 0.05;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eLOGS',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'ePOWER',iet_);
        elftv{it}{1} = 'K';
        elftp{it}{1} = 'B';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('T')
            ename = ['LOGC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eLOGS';
            ielftype(ie) = iet_('eLOGS');
            vname = ['C',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('C',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['KS',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePOWER';
            ielftype(ie) = iet_('ePOWER');
            vname = ['K',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('K',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('B',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('B');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('T')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['LOGC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['BT',int2str(I)]);
            ig = ig_(['NL',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['KS',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['AT',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -9.7457259D-01
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-V-V';
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

    case 'eLOGS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = log(EV_(1));
        if(nargout>1)
            g_(1,1) = 1.0/EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -1.0/EV_(1)^2;
                varargout{3} = H_;
            end
        end

    case 'ePOWER'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^pbm.elpar{iel_}(1);
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*EV_(1)^(pbm.elpar{iel_}(1)-1.0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) =...
                      pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0)*EV_(1)^(pbm.elpar{iel_}(1)-2.0);
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

