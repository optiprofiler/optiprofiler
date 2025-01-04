function varargout = HADAMALS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HADAMALS
%    --------
% 
%    An attempt to find Hadamard matrices of order N.
% 
%    The problem is to find an N by N orthonormal matrix Q,
%    with column norms N, whose entries are plus or minus one.
% 
%    Source:  A suggestion by Alan Edelman (MIT).
% 
%    SIF input: Nick Gould, Nov 1993.
% 
%    classification = 'C-COBR2-RN-V-V'
% 
%    The dimension of the matrix (=> N**2 variables).
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER    original value
% IE N                   4              $-PARAMETER
% IE N                   6              $-PARAMETER
% IE N                   8              $-PARAMETER
% IE N                   10             $-PARAMETER
% IE N                   12             $-PARAMETER
% IE N                   14             $-PARAMETER
% IE N                   16             $-PARAMETER
% IE N                   18             $-PARAMETER
% IE N                   20             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HADAMALS';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   32             $-PARAMETER
% IE N                   64             $-PARAMETER
% IE N                   128            $-PARAMETER
% IE N                   256            $-PARAMETER
% IE N                   428            $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        v_('N/2') = fix(v_('N')/v_('2'));
        v_('N/2+1') = 1+v_('N/2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['Q',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Q',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                [ig,ig_] = s2mpjlib('ii',['O',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        for J=v_('1'):v_('N')
            for I=v_('2'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['S',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for J=v_('1'):v_('N')
            pbm.gconst(ig_(['O',int2str(J),',',int2str(J)])) = v_('RN');
        end
        for J=v_('1'):v_('N')
            for I=v_('2'):v_('N')
                pbm.gconst(ig_(['S',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -1.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        for I=v_('1'):v_('N/2')
            pb.xlower(ix_(['Q',int2str(I),',',int2str(round(v_('1')))]),1) = 1.0;
            pb.xupper(ix_(['Q',int2str(I),',',int2str(round(v_('1')))]),1) = 1.0;
        end
        for I=v_('N/2+1'):v_('N')
            pb.xlower(ix_(['Q',int2str(I),',',int2str(round(v_('1')))]),1) = -1.0;
            pb.xupper(ix_(['Q',int2str(I),',',int2str(round(v_('1')))]),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N/2')
                pb.x0(ix_(['Q',int2str(I),',',int2str(J)]),1) = 0.9;
            end
            for I=v_('N/2+1'):v_('N')
                pb.x0(ix_(['Q',int2str(I),',',int2str(J)]),1) = -0.9;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'Q1';
        [it,iet_] = s2mpjlib( 'ii', 'en2PROD',iet_);
        elftv{it}{1} = 'Q1';
        elftv{it}{2} = 'Q2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                for K=v_('1'):v_('N')
                    ename = ['O',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PROD';
                    ielftype(ie) = iet_('en2PROD');
                    vname = ['Q',int2str(K),',',int2str(I)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1.0,1.0,[]);
                    posev = find(strcmp('Q1',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['Q',int2str(K),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1.0,1.0,[]);
                    posev = find(strcmp('Q2',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for J=v_('1'):v_('N')
            for I=v_('2'):v_('N')
                ename = ['S',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQR';
                ielftype(ie) = iet_('eSQR');
                vname = ['Q',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-1.0,1.0,[]);
                posev = find(strcmp('Q1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        [it,igt_] = s2mpjlib('ii','gLARGEL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('1'):v_('N')
            for I=v_('1'):J
                ig = ig_(['O',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gL2';
                for K=v_('1'):v_('N')
                    ig = ig_(['O',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['O',int2str(I),',',int2str(J),',',int2str(K)]);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        for J=v_('1'):v_('N')
            for I=v_('2'):v_('N')
                ig = ig_(['S',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gLARGEL2';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COBR2-RN-V-V';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e+0;
                varargout{3} = H_;
            end
        end

    case 'en2PROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0e+0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'g_globs'

        pbm = varargin{1};
        pbm.gfpar(1) = 1.0e+0;    % this is FACTOR
        varargout{1} = pbm;

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0e+0;
                varargout{3} = H_;
            end
        end

    case 'gLARGEL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = pbm.gfpar(1)*GVAR_*GVAR_;
        if(nargout>1)
            g_ = 2.0e+0*pbm.gfpar(1)*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0e+0*pbm.gfpar(1);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,1];
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

