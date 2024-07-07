function varargout = REPEAT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : REPEAT
%    *********
% 
%    This problem is to find the nearest feasible point to 2n+1 inconsistent
%    linear equations subject to bounds
% 
%    Source: blue-cheese delerium
% 
%    SIF input: Nick Gould, December 2020.
% 
%    classification = 'NLR2-AN-V-V'
% 
%    N is the number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'REPEAT';

switch(action)

    case 'setup'

    pb.name      = 'REPEAT';
    pb.sifpbname = 'REPEAT';
    pbm.name     = 'REPEAT';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('N-1') = -1+v_('N');
        v_('1') = 1;
        v_('2') = 2;
        v_('100') = 100;
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N/100') = fix(v_('N')/v_('100'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2xlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['C',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        [ig,ig_] = s2xlib('ii',['C',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['C',int2str(round(v_('N')))];
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N-1')
            v_('I+1') = 1+I;
            [ig,ig_] = s2xlib('ii',['R',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        [ig,ig_] = s2xlib('ii',['R',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['R',int2str(round(v_('N')))];
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N/100'):v_('N')
            v_('RI') = I;
            [ig,ig_] = s2xlib('ii','E',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'E';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('RI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('RI');
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
        pbm.congrps = find(ismember(gtype,{'<=','==','>='}));
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['C',int2str(I)])) = 2.0;
        end
        pbm.gconst(ig_(['C',int2str(round(v_('N')))])) = 1.0;
        for I=v_('1'):v_('N-1')
            pbm.gconst(ig_(['R',int2str(I)])) = 4.0;
        end
        pbm.gconst(ig_(['R',int2str(round(v_('N')))])) = 3.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = -1.0*ones(pb.n,1);
        pb.xupper(ix_(['X',int2str(round(v_('1')))])) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('2')))]),1) = 3.0;
        pb.xupper(ix_(['X',int2str(round(v_('N/2')))])) = 0.0;
        pb.xlower(ix_(['X',int2str(round(v_('N-1')))]),1) = 3.0;
        pb.xupper(ix_(['X',int2str(round(v_('N')))])) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'NLR2-AN-V-V';
        varargout{1} = pb;
        varargout{2} = pbm;

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
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

