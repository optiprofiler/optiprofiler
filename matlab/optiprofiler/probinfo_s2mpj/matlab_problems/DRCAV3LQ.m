function varargout = DRCAV3LQ(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DRCAV3LQ
%    *********
% 
%    This system of nonlinear equations models the stream function
%    corresponding to an incompressible fluid flow in a driven cavity 
%    (after elimination of the vorticity). The system is solved in the
%    least-squares sense.  The nonlinear system formulation is problem DRCAVTY3.
% 
%    The problem is nonconvex.
%    It differs from the problems DRCAV1LQ and DRCAV2LQ by the value 
%    chosen for the Reynolds number.
% 
%    Source:  
%    P.N. Brown and Y. Saad, 
%    "Hybrid Krylov Methods for Nonlinear Systems of Equations",
%    SIAM J. Sci. Stat. Comput. 11, pp. 450-481, 1990.
%    The boundary conditions have been set according to
%    I.E. Kaporin and O. Axelsson,
%    "On a class of nonlinear equation solvers based on the residual norm
%    reduction over a sequence of affine subspaces",
%    SIAM J, Sci. Comput. 16(1), 1995.
% 
%    SIF input: Ph. Toint, Jan 1995.
% 
%    classification = 'C-COXR2-MY-V-V'
% 
%    Discretization mesh: n = (M+3)**2 - fixed variables
% 
%       Alternative values for the SIF file parameters:
% IE M                   10             $-PARAMETER  n =   100   original value
% IE M                   31             $-PARAMETER  n =   961
% IE M                   63             $-PARAMETER  n =  3969
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DRCAV3LQ';

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
            v_('M') = 10;  %  SIF file default value
        else
            v_('M') = varargin{1};
        end
% IE M                   100            $-PARAMETER  n = 10000
        if(nargs<2)
            v_('RE') = 4500.0;  %  SIF file default value
        else
            v_('RE') = varargin{2};
        end
        v_('M+2') = 2+v_('M');
        v_('RM+2') = v_('M+2');
        v_('H') = 1.0/v_('RM+2');
        v_('-1') = -1;
        v_('0') = 0;
        v_('1') = 1;
        v_('M+1') = 1+v_('M');
        v_('H/2') = 0.5*v_('H');
        v_('-H/2') = -1.0*v_('H/2');
        v_('RE/4') = 0.25*v_('RE');
        v_('-RE/4') = -1.0*v_('RE/4');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('-1'):v_('M+2')
            for J=v_('-1'):v_('M+2')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Y',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            v_('I-2') = -2+I;
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            for J=v_('1'):v_('M')
                v_('J-2') = -2+J;
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                v_('J+2') = 2+J;
                [ig,ig_] = s2mpjlib('ii',['E',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
                iv = ix_(['Y',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 20.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 20.0;
                end
                iv = ix_(['Y',int2str(round(v_('I-1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -8.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -8.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -8.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -8.0;
                end
                iv = ix_(['Y',int2str(I),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -8.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -8.0;
                end
                iv = ix_(['Y',int2str(I),',',int2str(round(v_('J+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -8.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -8.0;
                end
                iv = ix_(['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 2.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 2.0;
                end
                iv = ix_(['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 2.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J+1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 2.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 2.0;
                end
                iv = ix_(['Y',int2str(round(v_('I-2'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+2'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(I),',',int2str(round(v_('J-2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(I),',',int2str(round(v_('J+2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        for J=v_('-1'):v_('M+2')
            pb.xlower(ix_(['Y',int2str(round(v_('-1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(round(v_('-1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['Y',int2str(round(v_('0'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(round(v_('0'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('1'):v_('M')
            pb.xlower(ix_(['Y',int2str(I),',',int2str(round(v_('-1')))]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(I),',',int2str(round(v_('-1')))]),1) = 0.0;
            pb.xlower(ix_(['Y',int2str(I),',',int2str(round(v_('0')))]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(I),',',int2str(round(v_('0')))]),1) = 0.0;
        end
        for I=v_('1'):v_('M')
            pb.xlower(ix_(['Y',int2str(I),',',int2str(round(v_('M+1')))]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(I),',',int2str(round(v_('M+1')))]),1) = 0.0;
            pb.xlower(ix_(['Y',int2str(I),',',int2str(round(v_('M+2')))]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(I),',',int2str(round(v_('M+2')))]),1) = 0.0;
        end
        for J=v_('-1'):v_('M+2')
            pb.xlower(ix_(['Y',int2str(round(v_('M+1'))),',',int2str(J)]),1) =...
                  v_('-H/2');
            pb.xupper(ix_(['Y',int2str(round(v_('M+1'))),',',int2str(J)]),1) =...
                  v_('-H/2');
            pb.xlower(ix_(['Y',int2str(round(v_('M+2'))),',',int2str(J)]),1) =...
                  v_('H/2');
            pb.xupper(ix_(['Y',int2str(round(v_('M+2'))),',',int2str(J)]),1) =...
                  v_('H/2');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eIPR',iet_);
        elftv{it}{1} = 'A1';
        elftv{it}{2} = 'A2';
        elftv{it}{3} = 'B1';
        elftv{it}{4} = 'B2';
        elftv{it}{5} = 'B3';
        elftv{it}{6} = 'B4';
        elftv{it}{7} = 'B5';
        elftv{it}{8} = 'B6';
        elftv{it}{9} = 'B7';
        elftv{it}{10} = 'B8';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('M')
            v_('I-2') = -2+I;
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            for J=v_('1'):v_('M')
                v_('J-2') = -2+J;
                v_('J-1') = -1+J;
                v_('J+1') = 1+J;
                v_('J+2') = 2+J;
                ename = ['X',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eIPR';
                    ielftype(ie) = iet_('eIPR');
                end
                vname = ['Y',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('A1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('A2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-2'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B5',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B6',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B7',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+2'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B8',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Z',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eIPR';
                    ielftype(ie) = iet_('eIPR');
                end
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('A1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('A2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J-2')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J-1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B4',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B5',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I-1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B6',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(round(v_('I+1'))),',',int2str(round(v_('J+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B7',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(I),',',int2str(round(v_('J+2')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('B8',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('M')
            for J=v_('1'):v_('M')
                ig = ig_(['E',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['X',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('RE/4');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Z',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_('-RE/4');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COXR2-MY-V-V';
        pb.x0          = zeros(pb.n,1);
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

    case 'eIPR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,10);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,3) = U_(2,3)+1;
        U_(2,4) = U_(2,4)+1;
        U_(2,5) = U_(2,5)+1;
        U_(2,6) = U_(2,6)-4;
        U_(2,7) = U_(2,7)+4;
        U_(2,8) = U_(2,8)-1;
        U_(2,9) = U_(2,9)-1;
        U_(2,10) = U_(2,10)-1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)*IV_(2);
        if(nargout>1)
            g_(1,1) = IV_(2);
            g_(2,1) = IV_(1);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = U_.'*H_*U_;
            end
        end

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

