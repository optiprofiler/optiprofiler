function varargout = TESTXSC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TESTXSC
%    *********
% 
%    A very simple quadratic problem designed to check the decoding of variable scaling.
%    The problem is to minimize
% 
%    f(x_1,...,x_n) = \sum_{i=1}^m [ \sum_{j=1,j odd} i x_j/j + \sum_{j=1,j
%     even} i x_j - 1 ]^2
% 
%    starting from x0 = (0, ..., 0).
% 
%    SIF input: Ph. Toint, Feb. 2024.
% 
%    classification = 'QUR2-AN-V-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TESTXSC';

switch(action)

    case 'setup'

    pb.name      = 'TESTXSC';
    pb.sifpbname = 'TESTXSC';
    pbm.name     = 'TESTXSC';
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        if(nargin<2)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        if(nargin<3)
            v_('M') = 1;  %  SIF file default value
        else
            v_('M') = varargin{2};
        end
%       Alternative values for the SIF file parameters:
% IE M                   N              $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('N')
            [iv,ix_] = s2xlib('ii',['X',int2str(J)],ix_);
            pb.xnames{iv} = ['X',int2str(J)];
        end
        for J=v_('1'):v_('2'):v_('N')
            v_('RJ') = J;
            [iv,ix_] = s2xlib('ii',['X',int2str(J)],ix_);
            pb.xnames{iv} = ['X',int2str(J)];
            xscale(iv,1) = v_('RJ');
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            v_('RI') = I;
            for J=v_('1'):v_('N')
                [ig,ig_] = s2xlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('RI')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('RI');
                end
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['G',int2str(I)])) = 1.0;
        end
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2xlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%%%% VARIABLES' SCALING %%%%%%%%%%%%%%%
        sA2 = size(pbm.A,2);
        for j = 1:min([sA2,pb.n,length(xscale)])
            if ( xscale(j) ~= 0.0 && xscale(j) ~= 1.0 )
                for i = find(pbm.A(:,j))
                      pbm.A(i,j) = pbm.A(i,j)/xscale(j);
                end
            end
        end
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'QUR2-AN-V-0';
        pb.x0          = zeros(pb.n,1);
        varargout{1} = pb;
        varargout{2} = pbm;

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

