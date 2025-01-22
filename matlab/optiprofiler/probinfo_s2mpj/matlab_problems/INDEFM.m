function varargout = INDEFM(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : INDEFM
%    *********
% 
%    Variant of INDEF, a nonconvex problem which has an indefinite Hessian 
%    at the starting point, by Luksan et al
% 
%    Source: problem 37 in
%    L. Luksan, C. Matonoha and J. Vlcek  
%    Modified CUTE problems for sparse unconstraoined optimization
%    Technical Report 1081
%    Institute of Computer Science
%    Academy of Science of the Czech Republic
% 
%    based on the original problem by N. Gould
% 
%    SIF input: Nick Gould, June, 2013
% 
%    classification = 'C-COUR2-AN-V-0'
% 
%    The number of variables is N.
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'INDEFM';

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
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   5000           $-PARAMETER
% IE N                   10000          $-PARAMETER
% IE N                   100000         $-PARAMETER
        if(nargs<2)
            v_('ALPHA') = 0.5;  %  SIF file default value
        else
            v_('ALPHA') = varargin{2};
        end
% RE ALPHA               1.0            $-PARAMETER
% RE ALPHA               10.0           $-PARAMETER
% RE ALPHA               100.0          $-PARAMETER
% RE ALPHA               1000.0         $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        v_('RN+1') = v_('N+1');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['SIN',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
        end
        for I=v_('2'):v_('N-1')
            [ig,ig_] = s2mpjlib('ii',['COS',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 2.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('N')))]);
            valA(end+1) = -1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('1')))]);
            valA(end+1) = -1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('T') = v_('RI')/v_('RN+1');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('T');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gCOS',igt_);
        [it,igt_] = s2mpjlib('ii','gCOS',igt_);
        grftp{it}{1} = 'ALPHA';
        [it,igt_] = s2mpjlib('ii','gSIN',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('2'):v_('N-1')
            ig = ig_(['COS',int2str(I)]);
            pbm.grftype{ig} = 'gCOS';
            [~,posgp] = ismember('ALPHA',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = v_('ALPHA');
        end
        for I=v_('1'):v_('N')
            ig = ig_(['SIN',int2str(I)]);
            pbm.grftype{ig} = 'gSIN';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               ??
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COUR2-AN-V-0';
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

% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gCOS'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = pbm.grpar{igr_}(1)*cos(GVAR_);
        if(nargout>1)
            g_ = -pbm.grpar{igr_}(1)*sin(GVAR_);
            varargout{2} = g_;
            if(nargout>2)
                H_ = -pbm.grpar{igr_}(1)*cos(GVAR_);
                varargout{3} = H_;
            end
        end

    case 'gSIN'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 100.0*sin(0.01*GVAR_);
        if(nargout>1)
            g_ = cos(0.01*GVAR_);
            varargout{2} = g_;
            if(nargout>2)
                H_ = -0.01*sin(0.01*GVAR_);
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

