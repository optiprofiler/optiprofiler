function varargout = CHNROSNB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHNROSNB
%    --------
%    The chained Rosenbrock function (Toint)
% 
%    Source:
%    Ph.L. Toint,
%    "Some numerical results using a sparse matrix updating formula in
%    unconstrained optimization",
%    Mathematics of Computation, vol. 32(114), pp. 839-852, 1978.
% 
%    See also Buckley#46 (n = 25) (p. 45).
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SUR2-AN-V-0'
% 
%    Number of variables ( at most 50)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   25             $-PARAMETER
% IE N                   50             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHNROSNB';

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
        v_('ALPH1') = 1.25;
        v_('ALPH2') = 1.40;
        v_('ALPH3') = 2.40;
        v_('ALPH4') = 1.40;
        v_('ALPH5') = 1.75;
        v_('ALPH6') = 1.20;
        v_('ALPH7') = 2.25;
        v_('ALPH8') = 1.20;
        v_('ALPH9') = 1.00;
        v_('ALPH10') = 1.10;
        v_('ALPH11') = 1.50;
        v_('ALPH12') = 1.60;
        v_('ALPH13') = 1.25;
        v_('ALPH14') = 1.25;
        v_('ALPH15') = 1.20;
        v_('ALPH16') = 1.20;
        v_('ALPH17') = 1.40;
        v_('ALPH18') = 0.50;
        v_('ALPH19') = 0.50;
        v_('ALPH20') = 1.25;
        v_('ALPH21') = 1.80;
        v_('ALPH22') = 0.75;
        v_('ALPH23') = 1.25;
        v_('ALPH24') = 1.40;
        v_('ALPH25') = 1.60;
        v_('ALPH26') = 2.00;
        v_('ALPH27') = 1.00;
        v_('ALPH28') = 1.60;
        v_('ALPH29') = 1.25;
        v_('ALPH30') = 2.75;
        v_('ALPH31') = 1.25;
        v_('ALPH32') = 1.25;
        v_('ALPH33') = 1.25;
        v_('ALPH34') = 3.00;
        v_('ALPH35') = 1.50;
        v_('ALPH36') = 2.00;
        v_('ALPH37') = 1.25;
        v_('ALPH38') = 1.40;
        v_('ALPH39') = 1.80;
        v_('ALPH40') = 1.50;
        v_('ALPH41') = 2.20;
        v_('ALPH42') = 1.40;
        v_('ALPH43') = 1.50;
        v_('ALPH44') = 1.25;
        v_('ALPH45') = 2.00;
        v_('ALPH46') = 1.50;
        v_('ALPH47') = 1.25;
        v_('ALPH48') = 1.40;
        v_('ALPH49') = 0.60;
        v_('ALPH50') = 1.50;
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('2'):v_('N')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['SQ',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            v_('AI2') = v_(['ALPH',int2str(I)])*v_(['ALPH',int2str(I)]);
            v_('16AI2') = 16.0*v_('AI2');
            v_('SCL') = 1.0/v_('16AI2');
            pbm.gscale(ig,1) = v_('SCL');
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('2'):v_('N')
            pbm.gconst(ig_(['B',int2str(I)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = -1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eETYPE',iet_);
        elftv{it}{1} = 'V1';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('2'):v_('N')
            ename = ['ELA',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eETYPE';
                ielftype(ie) = iet_('eETYPE');
            end
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],-1.0);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
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
        for I=v_('2'):v_('N')
            ig = ig_(['SQ',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ELA',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
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
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eETYPE'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = -EV_(1)^2;
        if(nargout>1)
            g_(1,1) = -2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -2.0;
                varargout{3} = H_;
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

