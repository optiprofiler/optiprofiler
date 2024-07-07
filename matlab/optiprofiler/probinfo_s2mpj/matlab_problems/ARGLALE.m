function varargout = ARGLALE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ARGLALE
%    *********
%    Variable dimension full rank linear problem, a linear equation version.
% 
%    Source: Problem 32 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#80 (with different N and M)
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'NLR2-AN-V-V'
% 
%    This is a(n infeasible) linear feasibility problem 
% 
%    N is the number of free variables
%    M is the number of equations ( M.ge.N)
% 
%       Alternative values for the SIF file parameters:
% IE N                   4              $-PARAMETER
% IE N                   10             $-PARAMETER
% IE N                   50             $-PARAMETER 
% IE N                   100            $-PARAMETER
% IE N                   200            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ARGLALE';

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
            v_('N') = 4;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE M                   6              $-PARAMETER .ge. N
% IE M                   20             $-PARAMETER .ge. N
% IE M                   100            $-PARAMETER .ge. N
% IE M                   200            $-PARAMETER .ge. N
% IE M                   400            $-PARAMETER .ge. N
        if(nargs<2)
            v_('M') = 6;  %  SIF file default value
        else
            v_('M') = varargin{2};
        end
        v_('1') = 1;
        v_('-2.0') = -2.0;
        v_('N+1') = 1+v_('N');
        v_('RM') = v_('M');
        v_('-2/M') = v_('-2.0')/v_('RM');
        v_('1-2/M') = 1.0+v_('-2/M');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            for J=v_('1'):v_('I-1')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/M')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/M');
                end
            end
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['G',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1-2/M')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1-2/M');
            end
            for J=v_('I+1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/M')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/M');
                end
            end
        end
        for I=v_('N+1'):v_('M')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/M')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/M');
                end
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
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = 1.0*ones(ngrp,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'NLR2-AN-V-V';
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

