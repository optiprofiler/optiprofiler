function varargout = KISSING(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem: KISSING NUMBER PROBLEM
%                                                                    
%    Source: This problem is associated to the family of Hard-Spheres 
%    problem. It belongs to the family of sphere packing problems, a 
%    class of challenging problems dating from the beginning of the 
%    17th century which is related to practical problems in Chemistry, 
%    Biology and Physics. It consists on maximizing the minimum pairwise 
%    distance between NP points on a sphere in \R^{MDIM}. 
%    This problem may be reduced to a nonconvex nonlinear optimization 
%    problem with a potentially large number of (nonoptimal) points 
%    satisfying optimality conditions. We have, thus, a class of problems 
%    indexed by the parameters MDIM and NP, that provides a suitable 
%    set of test problems for evaluating nonlinear programming codes.
%    After some algebric manipulations, we can formulate this problem as
%                             Minimize z
%                             subject to
%        
%       z \geq <x_i, x_j> for all different pair of indices i, j
%       
%                             ||x_i||^2 = 1    for all i = 1,...,NP
%      The goal is to find an objective value less than 0.5 (This means
%      that the NP points stored belong to the sphere and every distance
%      between two of them is greater than 1.0).
%      Obs: the starting point is aleatorally chosen although each 
%      variable belongs to [-1.,1.].
%      References:
%      [1] "Validation of an Augmented Lagrangian algorithm with a 
%           Gauss-Newton Hessian approximation using a set of 
%           Hard-Spheres problems", N. Krejic, J. M. Martinez, M. Mello 
%           and E. A. Pilotta, Tech. Report RP 29/98, IMECC-UNICAMP, 
%           Campinas, 1998.
%      [2] "Inexact-Restoration Algorithm for Constrained Optimization",
%           J. M. Martinez and E. A. Pilotta, Tech. Report, IMECC-UNICAMP, 
%           Campinas, 1998.
%      [3]  "Sphere Packings, Lattices and Groups", J. H. Conway and 
%            N. J. C. Sloane, Springer-Verlag, NY, 1988.
%      SIF input: September 29, 1998
% 		 Jose Mario Martinez
%                 Elvio Angel Pilotta
% 
%    classification = 'LQR2-RN-V-V'
% 
% **********************************************************************
% 
%    Number of points: NP >= 12
% 
%       Alternative values for the SIF file parameters:
% IE NP                   12            $-PARAMETER
% IE NP                   13            $-PARAMETER
% IE NP                   14            $-PARAMETER
% IE NP                   15            $-PARAMETER
% IE NP                   22            $-PARAMETER
% IE NP                   23            $-PARAMETER
% IE NP                   24            $-PARAMETER
% IE NP                   25            $-PARAMETER
% IE NP                   26            $-PARAMETER
% IE NP                   27            $-PARAMETER
% IE NP	                 37            $-PARAMETER
% IE NP                   38            $-PARAMETER
% IE NP                   39            $-PARAMETER
% IE NP                   40            $-PARAMETER
% IE NP                   41            $-PARAMETER
% IE NP                   42            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'KISSING';

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
            v_('NP') = 25;  %  SIF file default value
        else
            v_('NP') = varargin{1};
        end
% IE MDIM                 3             $-PARAMETER
        if(nargs<2)
            v_('MDIM') = 3;  %  SIF file default value
        else
            v_('MDIM') = varargin{2};
        end
% IE MDIM                 4             $-PARAMETER
% IE MDIM                 5             $-PARAMETER
        v_('N-') = -1+v_('NP');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('NP')
            for J=v_('1'):v_('MDIM')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        [iv,ix_] = s2mpjlib('ii','Z',ix_);
        pb.xnames{iv} = 'Z';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('Z');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('N-')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP')
                [ig,ig_] = s2mpjlib('ii',['IC',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['IC',int2str(I),',',int2str(J)];
                iv = ix_('Z');
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for I=v_('1'):v_('NP')
            [ig,ig_] = s2mpjlib('ii',['EC',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['EC',int2str(I)];
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('NP')
            pbm.gconst(ig_(['EC',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('NP')
            for J=v_('1'):v_('MDIM')
                pb.xlower(ix_(['X',int2str(I),',',int2str(J)])) = -Inf;
                pb.xupper(ix_(['X',int2str(I),',',int2str(J)]),1) = +Inf;
            end
        end
        pb.xlower(ix_('Z')) = -Inf;
        pb.xupper(ix_('Z'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('X1,1'),1) = -0.10890604;
        pb.x0(ix_('X1,2'),1) = 0.85395078;
        pb.x0(ix_('X1,3'),1) = -0.45461680;
        pb.x0(ix_('X2,1'),1) = 0.49883922;
        pb.x0(ix_('X2,2'),1) = -0.18439316;
        pb.x0(ix_('X2,3'),1) = -0.04798594;
        pb.x0(ix_('X3,1'),1) = 0.28262888;
        pb.x0(ix_('X3,2'),1) = -0.48054070;
        pb.x0(ix_('X3,3'),1) = 0.46715332;
        pb.x0(ix_('X4,1'),1) = -0.00580106;
        pb.x0(ix_('X4,2'),1) = -0.49987584;
        pb.x0(ix_('X4,3'),1) = -0.44130302;
        pb.x0(ix_('X5,1'),1) = 0.81712540;
        pb.x0(ix_('X5,2'),1) = -0.36874258;
        pb.x0(ix_('X5,3'),1) = -0.68321896;
        pb.x0(ix_('X6,1'),1) = 0.29642426;
        pb.x0(ix_('X6,2'),1) = 0.82315508;
        pb.x0(ix_('X6,3'),1) = 0.35938150;
        pb.x0(ix_('X7,1'),1) = 0.09215152;
        pb.x0(ix_('X7,2'),1) = -0.53564686;
        pb.x0(ix_('X7,3'),1) = 0.00191436;
        pb.x0(ix_('X8,1'),1) = 0.11700318;
        pb.x0(ix_('X8,2'),1) = 0.96722760;
        pb.x0(ix_('X8,3'),1) = -0.14916438;
        pb.x0(ix_('X9,1'),1) = 0.01791524;
        pb.x0(ix_('X9,2'),1) = 0.17759446;
        pb.x0(ix_('X9,3'),1) = -0.61875872;
        pb.x0(ix_('X10,1'),1) = -0.63833630;
        pb.x0(ix_('X10,2'),1) = 0.80830972;
        pb.x0(ix_('X10,3'),1) = 0.45846734;
        pb.x0(ix_('X11,1'),1) = 0.28446456;
        pb.x0(ix_('X11,2'),1) = 0.45686938;
        pb.x0(ix_('X11,3'),1) = 0.16368980;
        pb.x0(ix_('X12,1'),1) = 0.76557382;
        pb.x0(ix_('X12,2'),1) = 0.16700944;
        pb.x0(ix_('X12,3'),1) = -0.31647534;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eQUA',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N-')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP')
                for K=v_('1'):v_('MDIM')
                    ename = ['A',int2str(I),',',int2str(J),',',int2str(K)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'ePROD';
                    ielftype(ie) = iet_('ePROD');
                    vname = ['X',int2str(I),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('X',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(J),',',int2str(K)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('Y',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NP')
            for K=v_('1'):v_('MDIM')
                ename = ['B',int2str(I),',',int2str(K)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eQUA';
                ielftype(ie) = iet_('eQUA');
                vname = ['X',int2str(I),',',int2str(K)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N-')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP')
                for K=v_('1'):v_('MDIM')
                    ig = ig_(['IC',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J),',',int2str(K)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        for I=v_('1'):v_('NP')
            for K=v_('1'):v_('MDIM')
                ig = ig_(['EC',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(K)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION             4.47214D-01
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-RN-V-V';
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

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eQUA'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

