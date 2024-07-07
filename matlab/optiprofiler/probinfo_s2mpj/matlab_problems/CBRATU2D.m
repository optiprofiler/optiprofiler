function varargout = CBRATU2D(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CBRATU2D
%    *********
% 
%    The complex 2D Bratu problem on the unit square, using finite
%    differences.
% 
%    Source: problem 3 in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'NOR2-MN-V-V'
% 
%    P is the number of points in one side of the unit square (variable).
%    There are 2*P**2 variables
% 
%       Alternative values for the SIF file parameters:
% IE P                   4              $-PARAMETER n = 32     original value
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CBRATU2D';

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
% IE P                   7              $-PARAMETER n = 98
% IE P                   16             $-PARAMETER n = 512
% IE P                   23             $-PARAMETER n = 1058
% IE P                   40             $-PARAMETER n = 3200
        if(nargs<2)
            v_('LAMBDA') = 5.0;  %  SIF file default value
        else
            v_('LAMBDA') = varargin{2};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('1.0') = 1.0;
        v_('P-1') = -1+v_('P');
        v_('RP-1') = v_('P-1');
        v_('H') = v_('1.0')/v_('RP-1');
        v_('H2') = v_('H')*v_('H');
        v_('C') = v_('H2')*v_('LAMBDA');
        v_('-C') = -1.0*v_('C');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('P')
            for I=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('P-1')
            v_('R') = 1+I;
            v_('S') = -1+I;
            for J=v_('2'):v_('P-1')
                v_('V') = 1+J;
                v_('W') = -1+J;
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I),',',int2str(J)];
                iv = ix_(['U',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 4.0;
                end
                iv = ix_(['U',int2str(round(v_('R'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['U',int2str(round(v_('S'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['U',int2str(I),',',int2str(round(v_('V')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['U',int2str(I),',',int2str(round(v_('W')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2mpjlib('ii',['F',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['F',int2str(I),',',int2str(J)];
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 4.0;
                end
                iv = ix_(['X',int2str(round(v_('R'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('S'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('V')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('W')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
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
        for J=v_('1'):v_('P')
            pb.xlower(ix_(['U',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['U',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('1'))),',',int2str(J)]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(round(v_('P'))),',',int2str(J)]),1) = 0.0;
        end
        for I=v_('2'):v_('P-1')
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xlower(ix_(['U',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['U',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('P')))]),1) = 0.0;
            pb.xlower(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('1')))]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eRPART',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        [it,iet_] = s2mpjlib( 'ii', 'eCPART',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                ename = ['A',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eRPART';
                ielftype(ie) = iet_('eRPART');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['B',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eCPART';
                ielftype(ie) = iet_('eCPART');
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('V',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('2'):v_('P-1')
            for J=v_('2'):v_('P-1')
                ig = ig_(['G',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['A',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-C');
                ig = ig_(['F',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['B',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-C');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-V-V';
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

    case 'eRPART'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPU = exp(EV_(1));
        EXPUC = EXPU*cos(EV_(2));
        EXPUS = EXPU*sin(EV_(2));
        varargout{1} = EXPUC;
        if(nargout>1)
            g_(1,1) = EXPUC;
            g_(2,1) = -EXPUS;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = EXPUC;
                H_(1,2) = -EXPUS;
                H_(2,1) = H_(1,2);
                H_(2,2) = -EXPUC;
                varargout{3} = H_;
            end
        end

    case 'eCPART'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPU = exp(EV_(1));
        EXPUC = EXPU*cos(EV_(2));
        EXPUS = EXPU*sin(EV_(2));
        varargout{1} = EXPUS;
        if(nargout>1)
            g_(1,1) = EXPUS;
            g_(2,1) = EXPUC;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = EXPUS;
                H_(1,2) = EXPUC;
                H_(2,1) = H_(1,2);
                H_(2,2) = -EXPUS;
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

