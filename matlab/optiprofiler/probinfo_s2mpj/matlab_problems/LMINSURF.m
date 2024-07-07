function varargout = LMINSURF(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LMINSURF
%    *********
% 
%    The linear minimum surface problem.
% 
%    The problem comes from the discretization of the minimum surface
%    problem on the unit square: given a set of boundary conditions on
%    the four sides of the square, one must find the surface which
%    meets these boundary conditions and is of minimum area.
% 
%    The unit square is discretized into (p-1)**2 little squares. The
%    heights of the considered surface above the corners of these little
%    squares are the problem variables,  There are p**2 of them.
%    Given these heights, the area above a little square is
%    approximated by the
%      S(i,j) = sqrt( 1 + 0.5(p-1)**2 ( a(i,j)**2 + b(i,j)**2 ) ) / (p-1)**2
%    where
%      a(i,j) = x(i,j) - x(i+1,j+1)
%    and
%      b(i,j) = x(i+1,j) - x(i,j+1)
% 
%    In the Linear Mininum Surface, the boundary conditions are given
%    as the heights of a given plane above the square boundaries.  This
%    plane is specified by its height above the (0,0) point (H00 below),
%    and its slopes along the first and second coordinate
%    directions in the plane (these slopes are denoted SLOPEJ and SLOPEI below).
% 
%    Source:
%    A Griewank and Ph. Toint,
%    "Partitioned variable metric updates for large structured
%    optimization problems",
%    Numerische Mathematik 39:429-448, 1982.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'OXR2-MY-V-0'
% 
%    P is the number of points in one side of the unit square
% 
%       Alternative values for the SIF file parameters:
% IE P                   4              $-PARAMETER n = 16     original value
% IE P                   7              $-PARAMETER n = 49
% IE P                   8              $-PARAMETER n = 64
% IE P                   11             $-PARAMETER n = 121
% IE P                   31             $-PARAMETER n = 961
% IE P                   32             $-PARAMETER n = 1024
% IE P                   75             $-PARAMETER n = 5625
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LMINSURF';

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
            v_('P') = 4;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
% IE P                   100            $-PARAMETER n = 10000
% IE P                   125            $-PARAMETER n = 15625
        v_('H00') = 1.0;
        v_('SLOPEJ') = 4.0;
        v_('SLOPEI') = 8.0;
        v_('TWOP') = v_('P')+v_('P');
        v_('P-1') = -1+v_('P');
        v_('PP-1') = v_('P')*v_('P-1');
        v_('RP-1') = v_('P-1');
        v_('INVP-1') = 1.0/v_('RP-1');
        v_('RP-1SQ') = v_('INVP-1')*v_('INVP-1');
        v_('SCALE') = 1.0/v_('RP-1SQ');
        v_('SQP-1') = v_('RP-1')*v_('RP-1');
        v_('PARAM') = 0.5*v_('SQP-1');
        v_('1') = 1;
        v_('2') = 2;
        v_('STON') = v_('INVP-1')*v_('SLOPEI');
        v_('WTOE') = v_('INVP-1')*v_('SLOPEJ');
        v_('H01') = v_('H00')+v_('SLOPEJ');
        v_('H10') = v_('H00')+v_('SLOPEI');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('P')
            for I=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('P-1')
            for J=v_('1'):v_('P-1')
                [ig,ig_] = s2mpjlib('ii',['S',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                pbm.gscale(ig,1) = v_('SCALE');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
        pbm.gconst = -1.0*ones(ngrp,1);
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for J=v_('1'):v_('P')
            v_('J-1') = -1+J;
            v_('RJ-1') = v_('J-1');
            v_('TH') = v_('RJ-1')*v_('WTOE');
            v_('TL') = v_('TH')+v_('H00');
            v_('TU') = v_('TH')+v_('H10');
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = v_('TL');
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = v_('TL');
            pb.xlower(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = v_('TU');
            pb.xupper(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = v_('TU');
        end
        for I=v_('2'):v_('P-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TV') = v_('RI-1')*v_('STON');
            v_('TR') = v_('TV')+v_('H00');
            v_('TL') = v_('TV')+v_('H01');
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = v_('TL');
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = v_('TL');
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = v_('TR');
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = v_('TR');
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        for J=v_('1'):v_('P')
            v_('J-1') = -1+J;
            v_('RJ-1') = v_('J-1');
            v_('TH') = v_('RJ-1')*v_('WTOE');
            v_('TL') = v_('TH')+v_('H00');
            v_('TU') = v_('TH')+v_('H10');
            pb.x0(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = v_('TL');
            pb.x0(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = v_('TU');
        end
        for I=v_('2'):v_('P-1')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_('TV') = v_('RI-1')*v_('STON');
            v_('TR') = v_('TV')+v_('H00');
            v_('TL') = v_('TV')+v_('H01');
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = v_('TL');
            pb.x0(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = v_('TR');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eISQ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('P-1')
            v_('I+1') = 1+I;
            for J=v_('1'):v_('P-1')
                v_('J+1') = 1+J;
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eISQ';
                ielftype(ie) = iet_('eISQ');
                vname = ['X',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gSQROOT',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('P-1')
            for J=v_('1'):v_('P-1')
                ig = ig_(['S',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gSQROOT';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('PARAM');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('PARAM');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               9.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OXR2-MY-V-0';
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

    case 'eISQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        varargout{1} = IV_(1)*IV_(1);
        if(nargout>1)
            g_(1,1) = IV_(1)+IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gSQROOT'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        SQRAL = sqrt(GVAR_);
        varargout{1} = SQRAL;
        if(nargout>1)
            g_ = 0.5e0/SQRAL;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -0.25e0/(SQRAL*GVAR_);
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

