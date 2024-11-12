function varargout = CHENHARK(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHENHARK
%    --------
% 
%    A bound-constrained version the Linear Complementarity problem
% 
%    Find x such that w = M x + q, x and w nonnegative and x^T w = 0,
%    where
% 
%    M = (  6   -4   1   0  ........ 0 ) 
%        ( -4    6  -4   1  ........ 0 )
%        (  1   -4   6  -4  ........ 0 )
%        (  0    1  -4   6  ........ 0 )  
%           ..........................
%        (  0   ........... 0  1 -4  6 )
% 
%    and q is given.
% 
%    Source: 
%    B. Chen and P. T. Harker,
%    SIMAX 14 (1993) 1168-1190
% 
%    SDIF input: Nick Gould, November 1993.
% 
%    classification = 'C-CQBR2-AN-V-V'
% 
%    Number of variables
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   1000           $-PARAMETER     original value
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHENHARK';

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
        if(nargs<1)
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   10000          $-PARAMETER
% IE N                   50000          $-PARAMETER
% IE NFREE               5              $-PARAMETER
% IE NFREE               50             $-PARAMETER
% IE NFREE               500            $-PARAMETER     original value
% IE NFREE               2500           $-PARAMETER
        if(nargs<2)
            v_('NFREE') = 5;  %  SIF file default value
        else
            v_('NFREE') = varargin{2};
        end
% IE NFREE               5000           $-PARAMETER
% IE NFREE               10000          $-PARAMETER
% IE NDEGEN              2              $-PARAMETER
% IE NDEGEN              20             $-PARAMETER
% IE NDEGEN              200            $-PARAMETER     original value
% IE NDEGEN              500            $-PARAMETER
        if(nargs<3)
            v_('NDEGEN') = 2;  %  SIF file default value
        else
            v_('NDEGEN') = varargin{3};
        end
% IE NDEGEN              1000           $-PARAMETER
% IE NDEGEN              2000           $-PARAMETER
        v_('-1') = -1;
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('N-1') = -1+v_('N');
        v_('N+1') = 1+v_('N');
        v_('N+2') = 2+v_('N');
        v_('NFREE+1') = 1+v_('NFREE');
        v_('NF+ND') = v_('NFREE')+v_('NDEGEN');
        v_('NF+ND+1') = 1+v_('NF+ND');
        v_(['X',int2str(round(v_('-1')))]) = 0.0;
        v_(['X',int2str(round(v_('0')))]) = 0.0;
        for I=v_('1'):v_('NFREE')
            v_(['X',int2str(I)]) = 1.0;
        end
        for I=v_('NFREE+1'):v_('N+2')
            v_(['X',int2str(I)]) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('N-1')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['Q',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0;
            end
        end
        [ig,ig_] = s2mpjlib('ii',['Q',int2str(round(v_('0')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['Q',int2str(round(v_('1')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_(['X',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['Q',int2str(round(v_('N')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 2.0;
        end
        iv = ix_(['X',int2str(round(v_('N-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['Q',int2str(round(v_('N+1')))],ig_);
        gtype{ig} = '<>';
        iv = ix_(['X',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('NF+ND')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I-1') = -1+I;
            v_('I-2') = -2+I;
            v_('Q1') = -6.0*v_(['X',int2str(I)]);
            v_('Q2') = 4.0*v_(['X',int2str(round(v_('I+1')))]);
            v_('Q3') = 4.0*v_(['X',int2str(round(v_('I-1')))]);
            v_('Q4') = -1.0*v_(['X',int2str(round(v_('I+2')))]);
            v_('Q5') = -1.0*v_(['X',int2str(round(v_('I-2')))]);
            v_('Q') = v_('Q1')+v_('Q2');
            v_('Q') = v_('Q')+v_('Q3');
            v_('Q') = v_('Q')+v_('Q4');
            v_('Q') = v_('Q')+v_('Q5');
            [ig,ig_] = s2mpjlib('ii','L',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('Q')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('Q');
            end
        end
        for I=v_('NF+ND+1'):v_('N')
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I-1') = -1+I;
            v_('I-2') = -2+I;
            v_('Q1') = -6.0*v_(['X',int2str(I)]);
            v_('Q2') = 4.0*v_(['X',int2str(round(v_('I+1')))]);
            v_('Q3') = 4.0*v_(['X',int2str(round(v_('I-1')))]);
            v_('Q4') = -1.0*v_(['X',int2str(round(v_('I+2')))]);
            v_('Q5') = -1.0*v_(['X',int2str(round(v_('I-2')))]);
            v_('Q') = v_('Q1')+v_('Q2');
            v_('Q') = v_('Q')+v_('Q3');
            v_('Q') = v_('Q')+v_('Q4');
            v_('Q') = v_('Q')+v_('Q5');
            v_('Q') = 1.0+v_('Q');
            [ig,ig_] = s2mpjlib('ii','L',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('Q')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('Q');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gHALFL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('0'):v_('N+1')
            ig = ig_(['Q',int2str(I)]);
            pbm.grftype{ig} = 'gHALFL2';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 1.0;
%    Solution
% LO SOLTN               -0.5
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CQBR2-AN-V-V';
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

    case 'gHALFL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 5.0e-1*GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 1.0e+0;
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

