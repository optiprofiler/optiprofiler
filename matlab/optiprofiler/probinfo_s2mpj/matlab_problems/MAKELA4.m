function varargout = MAKELA4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MAKELA4
%    *********
% 
%    A nonlinear minmax problem in twenty variables.
% 
%    Source: 
%    M.M. Makela,
%    "Nonsmooth optimization",
%    Ph.D. thesis, Jyvaskyla University, 1990
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'LLR2-AN-21-40'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MAKELA4';

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
        v_('20') = 20;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('20')
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
        for I=v_('1'):v_('20')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['F',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['MF',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['MF',int2str(I)];
            iv = ix_('U');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 1.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 2.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 3.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 3.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 4.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 4.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 5.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 5.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 6.0;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 6.0;
        end
        if(isKey(ix_,'X7'))
            pb.x0(ix_('X7'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('X7')),1) = 7.0;
        end
        if(isKey(ix_,'X8'))
            pb.x0(ix_('X8'),1) = 8.0;
        else
            pb.y0(find(pbm.congrps==ig_('X8')),1) = 8.0;
        end
        if(isKey(ix_,'X9'))
            pb.x0(ix_('X9'),1) = 9.0;
        else
            pb.y0(find(pbm.congrps==ig_('X9')),1) = 9.0;
        end
        if(isKey(ix_,'X10'))
            pb.x0(ix_('X10'),1) = 10.0;
        else
            pb.y0(find(pbm.congrps==ig_('X10')),1) = 10.0;
        end
        if(isKey(ix_,'X11'))
            pb.x0(ix_('X11'),1) = -11.0;
        else
            pb.y0(find(pbm.congrps==ig_('X11')),1) = -11.0;
        end
        if(isKey(ix_,'X12'))
            pb.x0(ix_('X12'),1) = -12.0;
        else
            pb.y0(find(pbm.congrps==ig_('X12')),1) = -12.0;
        end
        if(isKey(ix_,'X13'))
            pb.x0(ix_('X13'),1) = -13.0;
        else
            pb.y0(find(pbm.congrps==ig_('X13')),1) = -13.0;
        end
        if(isKey(ix_,'X14'))
            pb.x0(ix_('X14'),1) = -14.0;
        else
            pb.y0(find(pbm.congrps==ig_('X14')),1) = -14.0;
        end
        if(isKey(ix_,'X15'))
            pb.x0(ix_('X15'),1) = -15.0;
        else
            pb.y0(find(pbm.congrps==ig_('X15')),1) = -15.0;
        end
        if(isKey(ix_,'X16'))
            pb.x0(ix_('X16'),1) = -16.0;
        else
            pb.y0(find(pbm.congrps==ig_('X16')),1) = -16.0;
        end
        if(isKey(ix_,'X17'))
            pb.x0(ix_('X17'),1) = -17.0;
        else
            pb.y0(find(pbm.congrps==ig_('X17')),1) = -17.0;
        end
        if(isKey(ix_,'X18'))
            pb.x0(ix_('X18'),1) = -18.0;
        else
            pb.y0(find(pbm.congrps==ig_('X18')),1) = -18.0;
        end
        if(isKey(ix_,'X19'))
            pb.x0(ix_('X19'),1) = -19.0;
        else
            pb.y0(find(pbm.congrps==ig_('X19')),1) = -19.0;
        end
        if(isKey(ix_,'X20'))
            pb.x0(ix_('X20'),1) = -20.0;
        else
            pb.y0(find(pbm.congrps==ig_('X20')),1) = -20.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-AN-21-40';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

