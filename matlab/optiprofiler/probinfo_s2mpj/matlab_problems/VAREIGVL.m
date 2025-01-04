function varargout = VAREIGVL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : VAREIGVL
%    *********
% 
%    The variational eigenvalue by Auchmuty.
%    This problems features a banded matrix of bandwidth 2M+1 = 9.
% 
%    This problem has N least-squares groups, each having a linear part
%    only and N nonlinear elements,
%    plus a least q-th power group having N nonlinear elements.
% 
%    Source: problem 1 in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
%               correction by Ph. Shott, January, 1995
%               and Nick Gould, December, 2019, May 2024
% 
%    classification = 'C-COUR2-AN-V-0'
% 
%    Number of variables -1 (variable)
% 
%       Alternative values for the SIF file parameters:
% IE N                   19             $-PARAMETER
% IE N                   49             $-PARAMETER     original value
% IE N                   99             $-PARAMETER
% IE N                   499            $-PARAMETER
% IE N                   999            $-PARAMETER
% IE N                   4999           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'VAREIGVL';

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
            v_('N') = 19;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE M                   4              $-PARAMETER  .le. N
% IE M                   5              $-PARAMETER  .le. N
        if(nargs<2)
            v_('M') = 6;  %  SIF file default value
        else
            v_('M') = varargin{2};
        end
        if(nargs<3)
            v_('Q') = 1.5;  %  SIF file default value
        else
            v_('Q') = varargin{3};
        end
        v_('1') = 1;
        v_('-1.0') = -1.0;
        v_('N+1') = 1+v_('N');
        v_('-M') = -1*v_('M');
        v_('M+1') = 1+v_('M');
        v_('N-M') = v_('N')+v_('-M');
        v_('N-M+1') = 1+v_('N-M');
        v_('N2') = v_('N')*v_('N');
        v_('RN2') = v_('N2');
        v_('-1/N2') = v_('-1.0')/v_('RN2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','MU',ix_);
        pb.xnames{iv} = 'MU';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('-I') = -1.0*v_('RI');
            v_('I+M') = I+v_('M');
            for J=v_('1'):v_('I+M')
                v_('RJ') = J;
                v_('IJ') = v_('RI')*v_('RJ');
                v_('SIJ') = sin(v_('IJ'));
                v_('J-I') = v_('RJ')+v_('-I');
                v_('J-ISQ') = v_('J-I')*v_('J-I');
                v_('ARG') = v_('J-ISQ')*v_('-1/N2');
                v_('EX') = exp(v_('ARG'));
                v_('AIJ') = v_('SIJ')*v_('EX');
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('AIJ');
            end
        end
        for I=v_('M+1'):v_('N-M')
            v_('RI') = I;
            v_('-I') = -1.0*v_('RI');
            v_('I-M') = I+v_('-M');
            v_('I+M') = I+v_('M');
            for J=v_('I-M'):v_('I+M')
                v_('RJ') = J;
                v_('IJ') = v_('RI')*v_('RJ');
                v_('SIJ') = sin(v_('IJ'));
                v_('J-I') = v_('RJ')+v_('-I');
                v_('J-ISQ') = v_('J-I')*v_('J-I');
                v_('ARG') = v_('J-ISQ')*v_('-1/N2');
                v_('EX') = exp(v_('ARG'));
                v_('AIJ') = v_('SIJ')*v_('EX');
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('AIJ');
            end
        end
        for I=v_('N-M+1'):v_('N')
            v_('RI') = I;
            v_('-I') = -1.0*v_('RI');
            v_('I-M') = I+v_('-M');
            for J=v_('I-M'):v_('N')
                v_('RJ') = J;
                v_('IJ') = v_('RI')*v_('RJ');
                v_('SIJ') = sin(v_('IJ'));
                v_('J-I') = v_('RJ')+v_('-I');
                v_('J-ISQ') = v_('J-I')*v_('J-I');
                v_('ARG') = v_('J-ISQ')*v_('-1/N2');
                v_('EX') = exp(v_('ARG'));
                v_('AIJ') = v_('SIJ')*v_('EX');
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig} = '<>';
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(J)]);
                valA(end+1) = v_('AIJ');
            end
        end
        [ig,ig_] = s2mpjlib('ii',['G',int2str(round(v_('N+1')))],ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.x0(ix_('MU'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'M';
        elftv{it}{2} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['P',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = 'MU';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('M',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['S',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQ';
            ielftype(ie) = iet_('eSQ');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gLQ',igt_);
        [it,igt_] = s2mpjlib('ii','gLQ',igt_);
        grftp{it}{1} = 'POWER';
        [it,igt_] = s2mpjlib('ii','gLQ2',igt_);
        [it,igt_] = s2mpjlib('ii','gLQ2',igt_);
        grftp{it}{1} = 'POWER';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['G',int2str(I)]);
            pbm.grftype{ig} = 'gLQ';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
            [~,posgp] = ismember('POWER',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = 2.0;
        end
        ig = ig_(['G',int2str(round(v_('N+1')))]);
        pbm.grftype{ig} = 'gLQ2';
        for I=v_('1'):v_('N')
            ig = ig_(['G',int2str(round(v_('N+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['S',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_(['G',int2str(round(v_('N+1')))]);
        [~,posgp] = ismember('POWER',grftp{igt_(pbm.grftype{ig})});
        pbm.grpar{ig}(posgp) = v_('Q');
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
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

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PR'

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

    case 'eSQ'

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

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gLQ'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        IPOWER = pbm.grpar{igr_}(1);
        PM1 = IPOWER-1;
        varargout{1} = GVAR_^IPOWER/pbm.grpar{igr_}(1);
        if(nargout>1)
            g_ = GVAR_^PM1;
            varargout{2} = g_;
            if(nargout>2)
                H_ = PM1*GVAR_^(IPOWER-2);
                varargout{3} = H_;
            end
        end

    case 'gLQ2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^pbm.grpar{igr_}(1)/pbm.grpar{igr_}(1);
        if(nargout>1)
            g_ = GVAR_^(pbm.grpar{igr_}(1)-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_ = (pbm.grpar{igr_}(1)-1.0e0)*GVAR_^(pbm.grpar{igr_}(1)-2.0e0);
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

