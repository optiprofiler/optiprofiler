function varargout = UBH5(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : UBH5
%    *********
% 
%    The problem is to minimize the integral of the control magnitude needed
%    to bring a vehicle, from given position and velocity, to the origin with
%    zero velocity in a fixed amount of time.  The controls are the components
%    of the vehicle acceleration. The discretization uses the trapezoidal rule.
%    This version of the problem is a variant of UBH1, where the cumulative
%    value of the objective is maintained as an additional state variable.
% 
%    The problem is convex.
% 
%    Source: unscaled problem 5 
%    (ODE = 1, CLS = 2, GRD = 1, MET = T, SEED = 0.) in
%    J.T. Betts and W.P. Huffman,
%    "Sparse Nonlinear Programming Test Problems (Release 1.0)",
%    Boeing Computer services, Seattle, July 1993.
% 
%    SIF input: Ph.L. Toint, October 1993.
% 
%    classification = 'LQR2-MN-V-V'
% 
%    Number of grid points
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER n=100, m=70    original value
% IE N                   100            $-PARAMETER n=1000, m=700
% IE N                   500            $-PARAMETER n=5000, m=3500
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'UBH5';

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
% IE N                   1000           $-PARAMETER n=10000, m=7000
% IE N                   2000           $-PARAMETER n=20000, m=14000
        v_('T0') = 0.0;
        v_('TF') = 1000.0;
        v_('RN') = v_('N');
        v_('TTIME') = v_('TF')-v_('T0');
        v_('K') = v_('TTIME')/v_('RN');
        v_('-K/2') = -0.5*v_('K');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('7')
            for T=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(I),',',int2str(T)],ix_);
                pb.xnames{iv} = ['Y',int2str(I),',',int2str(T)];
            end
        end
        for I=v_('1'):v_('3')
            for T=v_('0'):v_('N')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(T)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(T)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_(['Y',int2str(round(v_('7'))),',',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        for I=v_('1'):v_('3')
            v_('I+3') = 3+I;
            for T=v_('1'):v_('N')
                v_('T-1') = -1+T;
                [ig,ig_] = s2mpjlib('ii',['S',int2str(I),',',int2str(T)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['S',int2str(I),',',int2str(T)];
                iv = ix_(['Y',int2str(I),',',int2str(T)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(I),',',int2str(round(v_('T-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+3'))),',',int2str(round(v_('T-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-K/2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-K/2');
                end
                iv = ix_(['Y',int2str(round(v_('I+3'))),',',int2str(T)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-K/2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-K/2');
                end
            end
        end
        for I=v_('1'):v_('3')
            v_('I+3') = 3+I;
            for T=v_('1'):v_('N')
                v_('T-1') = -1+T;
                [ig,ig_] =...
                      s2mpjlib('ii',['S',int2str(round(v_('I+3'))),',',int2str(T)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['S',int2str(round(v_('I+3'))),',',int2str(T)];
                iv = ix_(['Y',int2str(round(v_('I+3'))),',',int2str(T)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['Y',int2str(round(v_('I+3'))),',',int2str(round(v_('T-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] =...
                      s2mpjlib('ii',['S',int2str(round(v_('I+3'))),',',int2str(T)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['S',int2str(round(v_('I+3'))),',',int2str(T)];
                iv = ix_(['U',int2str(I),',',int2str(round(v_('T-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-K/2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-K/2');
                end
                [ig,ig_] =...
                      s2mpjlib('ii',['S',int2str(round(v_('I+3'))),',',int2str(T)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['S',int2str(round(v_('I+3'))),',',int2str(T)];
                iv = ix_(['U',int2str(I),',',int2str(T)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-K/2')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-K/2');
                end
            end
        end
        for T=v_('1'):v_('N')
            v_('T-1') = -1+T;
            [ig,ig_] = s2mpjlib('ii',['S',int2str(round(v_('7'))),',',int2str(T)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['S',int2str(round(v_('7'))),',',int2str(T)];
            iv = ix_(['Y',int2str(round(v_('7'))),',',int2str(T)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_(['Y',int2str(round(v_('7'))),',',int2str(round(v_('T-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for I=v_('1'):v_('3')
            for T=v_('0'):v_('N')
                pb.xlower(ix_(['U',int2str(I),',',int2str(T)]),1) = -1.0;
                pb.xupper(ix_(['U',int2str(I),',',int2str(T)])) = 1.0;
            end
        end
        pb.xlower(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xupper(ix_(['Y',int2str(round(v_('1'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xlower(ix_(['Y',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xupper(ix_(['Y',int2str(round(v_('2'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xlower(ix_(['Y',int2str(round(v_('3'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xupper(ix_(['Y',int2str(round(v_('3'))),',',int2str(round(v_('0')))]),1) = 1000.0;
        pb.xlower(ix_(['Y',int2str(round(v_('4'))),',',int2str(round(v_('0')))]),1) = -...
             10.0;
        pb.xupper(ix_(['Y',int2str(round(v_('4'))),',',int2str(round(v_('0')))]),1) = -...
             10.0;
        pb.xlower(ix_(['Y',int2str(round(v_('5'))),',',int2str(round(v_('0')))]),1) = 10.0;
        pb.xupper(ix_(['Y',int2str(round(v_('5'))),',',int2str(round(v_('0')))]),1) = 10.0;
        pb.xlower(ix_(['Y',int2str(round(v_('6'))),',',int2str(round(v_('0')))]),1) = -...
             10.0;
        pb.xupper(ix_(['Y',int2str(round(v_('6'))),',',int2str(round(v_('0')))]),1) = -...
             10.0;
        pb.xlower(ix_(['Y',int2str(round(v_('7'))),',',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('7'))),',',int2str(round(v_('0')))]),1) = 0.0;
        for I=v_('1'):v_('6')
            pb.xlower(ix_(['Y',int2str(I),',',int2str(round(v_('N')))]),1) = 0.0;
            pb.xupper(ix_(['Y',int2str(I),',',int2str(round(v_('N')))]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for T=v_('0'):v_('N')
            for I=v_('1'):v_('3')
                ename = ['E',int2str(I),',',int2str(T)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = ['U',int2str(I),',',int2str(T)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for T=v_('1'):v_('N')
            v_('T-1') = -1+T;
            for I=v_('1'):v_('3')
                ig = ig_(['S',int2str(round(v_('7'))),',',int2str(T)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(round(v_('T-1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-K/2');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-K/2');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN(10)           1.14735202967
% LO SOLTN(100)          1.11631518169
% LO SOLTN(1000)         1.11598643493
% LO SOLTN(2000)         1.11587382445
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'LQR2-MN-V-V';
        pb.x0          = zeros(pb.n,1);
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

