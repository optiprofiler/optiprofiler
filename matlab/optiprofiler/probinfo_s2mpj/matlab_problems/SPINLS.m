function varargout = SPINLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SPINLS
%    *********
% 
%    Given n particles z_j = x_j + i * y_j in the complex plane,
%    determine their positions so that the equations
% 
%      z'_j = lambda z_j,
% 
%    where z_j = sum_k \j i / conj( z_j - z_k ) and i = sqrt(-1)
%    for some lamda = mu + i * omega
% 
%    A problem posed by Nick Trefethen
% 
%    SIF input: Nick Gould, June 2009
% 
%    classification = 'SUR2-AN-V-0'
%    Least-squares version of SPIN.SIF, Nick Gould, Jan 2020.
% 
%    Number of particles n
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER matrix dimension
% IE N                   3              $-PARAMETER matrix dimension
% IE N                   5              $-PARAMETER matrix dimension
% IE N                   10             $-PARAMETER matrix dimension
% IE N                   15             $-PARAMETER matrix dimension
% IE N                   20             $-PARAMETER matrix dimension
% IE N                   25             $-PARAMETER matrix dimension
% IE N                   30             $-PARAMETER matrix dimension
% IE N                   35             $-PARAMETER matrix dimension
% IE N                   40             $-PARAMETER matrix dimension
% IE N                   45             $-PARAMETER matrix dimension
% IE N                   50             $-PARAMETER matrix dimension
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SPINLS';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('RN') = v_('N');
        v_('PI/4') = atan(1.0);
        v_('2PI') = 8.0*v_('PI/4');
        v_('2PI/N') = v_('2PI')/v_('RN');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','MU',ix_);
        pb.xnames{iv} = 'MU';
        [iv,ix_] = s2mpjlib('ii','OMEGA',ix_);
        pb.xnames{iv} = 'OMEGA';
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['R',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['I',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                [ig,ig_] = s2mpjlib('ii',['M',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('2PII/N') = v_('2PI/N')*v_('RI');
            v_('C') = cos(v_('2PII/N'));
            v_('S') = sin(v_('2PII/N'));
            pb.x0(ix_(['X',int2str(I)]),1) = v_('C');
            pb.x0(ix_(['Y',int2str(I)]),1) = v_('S');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eRATIO',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        elftv{it}{3} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eSQR2',iet_);
        elftv{it}{1} = 'Y';
        elftv{it}{2} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['MX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'MU';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['MY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'MU';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['OX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'OMEGA';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['OY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'OMEGA';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                ename = ['RX',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eRATIO';
                ielftype(ie) = iet_('eRATIO');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['RY',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eRATIO';
                ielftype(ie) = iet_('eRATIO');
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                ename = ['V',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQR';
                ielftype(ie) = iet_('eSQR');
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['X',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQR2';
                ielftype(ie) = iet_('eSQR2');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Y',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQR2';
                ielftype(ie) = iet_('eSQR2');
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
                posev = find(strcmp('Z',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            ig = ig_(['R',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['MX',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['OY',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
            for J=v_('1'):v_('I-1')
                ig = ig_(['R',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RY',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.0;
            end
            for J=v_('I+1'):v_('N')
                ig = ig_(['R',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RY',int2str(J),',',int2str(I)]);
                pbm.grelw{ig}(posel) = -1.0;
            end
            ig = ig_(['I',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['MY',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['OX',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
            for J=v_('1'):v_('I-1')
                ig = ig_(['I',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RX',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = -1.0;
            end
            for J=v_('I+1'):v_('N')
                ig = ig_(['I',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['RX',int2str(J),',',int2str(I)]);
                pbm.grelw{ig}(posel) = 1.0;
            end
        end
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('I-1')
                ig = ig_(['M',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['V',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = -1.0;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['X',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.0;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Y',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.0;
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
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
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eSQR2'

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

    case 'eRATIO'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        U_(2,3) = U_(2,3)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        varargout{1} = IV_(1)/IV_(2)^2;
        if(nargout>1)
            g_(1,1) = 1.0/IV_(2)^2;
            g_(2,1) = -2.0*IV_(1)/IV_(2)^3;
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -2.0/IV_(2)^3;
                H_(2,1) = H_(1,2);
                H_(2,2) = 6.0*IV_(1)/IV_(2)^4;
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

