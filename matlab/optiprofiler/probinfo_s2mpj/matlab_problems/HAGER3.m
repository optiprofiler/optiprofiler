function varargout = HAGER3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HAGER3
%    *********
% 
%    A nonlinear optimal control problem, by W. Hager.
% 
%    Source: problem P3 in
%    W.W. Hager,
%    "Multiplier Methods for Nonlinear Optimal Control",
%    SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.
% 
%    SIF input: Ph. Toint, March 1991.
% 
%    classification = 'SLR2-AY-V-V'
% 
%    Number of discretized points in [0,1]
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   50             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   2500           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HAGER3';

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
            v_('N') = 5;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   5000           $-PARAMETER
        v_('N-1') = -1+v_('N');
        v_('RN') = v_('N');
        v_('H') = 1.0/v_('RN');
        v_('4/H') = 4.0*v_('RN');
        v_('8/H') = 8.0*v_('RN');
        v_('1/H-1/4') = -0.25+v_('RN');
        v_('-1/H') = -1.0*v_('RN');
        v_('-1/H-1/4') = -0.25+v_('-1/H');
        v_('0') = 0;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_('8/H');
            [ig,ig_] = s2mpjlib('ii',['U',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            pbm.gscale(ig,1) = v_('4/H');
        end
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['S',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['S',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H-1/4')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H-1/4');
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H-1/4')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H-1/4');
            end
            iv = ix_(['U',int2str(I)]);
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
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 1.0;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('0')))]),1) = 1.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eLINSQ',iet_);
        elftv{it}{1} = 'XA';
        elftv{it}{2} = 'XB';
        [it,iet_] = s2mpjlib( 'ii', 'eLINU',iet_);
        elftv{it}{1} = 'XA';
        elftv{it}{2} = 'XB';
        elftv{it}{3} = 'U';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            v_('I-1') = -1+I;
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eLINSQ';
            ielftype(ie) = iet_('eLINSQ');
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['UX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eLINU';
            ielftype(ie) = iet_('eLINU');
            vname = ['X',int2str(round(v_('I-1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('XB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_(['OX',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 0.625;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['UX',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['U',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           0.14103316577
% LO SOLTN(50)           0.14096413824
% LO SOLTN(100)          0.14096197328
% LO SOLTN(500)          0.14096128049
% LO SOLTN(1000)         0.14096125887
% LO SOLTN(5000)         0.14096125191
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'SLR2-AY-V-V';
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

    case 'eLINSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)+EV_(1)*EV_(2)+EV_(2)*EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)+EV_(2);
            g_(2,1) = 2.0*EV_(2)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0;
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eLINU'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(2,3);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        U_(2,3) = U_(2,3)+1;
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

