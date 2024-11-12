function varargout = LISWET4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LISWET4
%    *********
% 
%    A k-convex approximation problem posed as a 
%    convex quadratic problem, with variable dimensions.
% 
%    Formulation:
%    -----------
% 
%                 n+k             2
%    minimize 1/2 sum ( x  - c  )
%                 i=1    i    i
% 
%    subject to
% 
%                  k              k-i
%                 sum ( k ) ( -1 )    x     > 0
%                 i=0 ( i )            j+i  = 
% 
%    where c  = g( t ) + small perturbation, t  = (i-1)/(n+k-1)
%           i       i                         i 
% 
%    Case 4: g(t) = t**3
% 
%    NB. Perturbations are not random as Li and Swetits's 
%        random number generator is undefined.
% 
%    Source:
%    W. Li and J. Swetits,
%    "A Newton method for convex regression, data smoothing and
%    quadratic programming with bounded constraints",
%    SIAM J. Optimization 3 (3) pp 466-488, 1993.
% 
%    SIF input: Nick Gould, August 1994.
% 
%    classification = 'C-CQLR2-AN-V-V'
% 
%       Alternative values for the SIF file parameters:
% IE N                   100            $-PARAMETER 103 variables original value 
% IE K                   3              $-PARAMETER original value
% 
% IE N                   100            $-PARAMETER 104 variables    
% IE K                   4              $-PARAMETER
% 
% IE N                   100            $-PARAMETER 105 variables    
% IE K                   5              $-PARAMETER
% 
% IE N                   100            $-PARAMETER 106 variables    
% IE K                   6              $-PARAMETER
% 
% IE N                   400            $-PARAMETER 402 variables    
% IE K                   2              $-PARAMETER
% 
% IE N                   400            $-PARAMETER 403 variables    
% IE K                   3              $-PARAMETER
% 
% IE N                   2000           $-PARAMETER 2001 variables    
% IE K                   1              $-PARAMETER
% 
% IE N                   2000           $-PARAMETER 2002 variables    
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LISWET4';

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
            v_('N') = 50;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE K                   2              $-PARAMETER
        if(nargs<2)
            v_('K') = 3;  %  SIF file default value
        else
            v_('K') = varargin{2};
        end
% IE N                   10000          $-PARAMETER 10001 variables    
% IE K                   1              $-PARAMETER
% IE N                   10000          $-PARAMETER 10002 variables    
% IE K                   2              $-PARAMETER
        v_('0') = 0;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('HALF') = 0.5;
        v_('N+K') = v_('N')+v_('K');
        v_('N+K-1') = -1+v_('N+K');
        v_('RN+K-1') = v_('N+K-1');
        v_('CONST') = 0.0;
        v_(['B',int2str(round(v_('0')))]) = v_('ONE');
        for I=v_('1'):v_('K')
            v_('I-1') = -1+I;
            v_('RI') = I;
            v_(['B',int2str(I)]) = v_(['B',int2str(round(v_('I-1')))])*v_('RI');
        end
        v_(['C',int2str(round(v_('0')))]) = v_('ONE');
        v_('PLUSMINUS') = v_('ONE');
        for I=v_('1'):v_('K')
            v_('K-I') = v_('K')-I;
            v_('PLUSMINUS') = -1.0*v_('PLUSMINUS');
            v_(['C',int2str(I)]) =...
                  v_(['B',int2str(round(v_('K')))])/v_(['B',int2str(I)]);
            v_(['C',int2str(I)]) =...
                  v_(['C',int2str(I)])/v_(['B',int2str(round(v_('K-I')))]);
            v_(['C',int2str(I)]) = v_(['C',int2str(I)])*v_('PLUSMINUS');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N+K')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N+K')
            v_('I-1') = -1+I;
            v_('RI') = I;
            v_('RI-1') = v_('I-1');
            v_('TI') = v_('RI-1')/v_('RN+K-1');
            v_('TI**2') = v_('TI')*v_('TI');
            v_('GT') = v_('TI**2')*v_('TI');
            v_('RANDOM') = sin(v_('RI'));
            v_('RANDOM') = 0.1*v_('RANDOM');
            v_('CI') = v_('GT')+v_('RANDOM');
            v_('-CI') = -1.0*v_('CI');
            v_('-CI*CI') = v_('-CI')*v_('CI');
            v_('CONST') = v_('CONST')+v_('-CI*CI');
            [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-CI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-CI');
            end
        end
        for J=v_('1'):v_('N')
            v_('J+K') = J+v_('K');
            for I=v_('0'):v_('K')
                v_('J+K-I') = v_('J+K')-I;
                [ig,ig_] = s2mpjlib('ii',['CON',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['CON',int2str(J)];
                iv = ix_(['X',int2str(round(v_('J+K-I')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['C',int2str(I)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['C',int2str(I)]);
                end
            end
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
        v_('CONST') = v_('HALF')*v_('CONST');
        pbm.gconst(ig_('OBJ')) = v_('CONST');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N+K')
            ename = ['XSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N+K')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['XSQ',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQLR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
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

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = 5.0e-1*EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 1.0e+0;
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

