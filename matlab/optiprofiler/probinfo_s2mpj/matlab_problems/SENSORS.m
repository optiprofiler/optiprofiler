function varargout = SENSORS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SENSORS
%    *********
% 
%    A problem arising from two-dimensional optimal sensor placement
% 
%    Source:
%    H. Zhang and X. Wang,
%    "Optimal sensor placement",
%    SIAM Review, vol. 35, p. 641, 1993.
% 
%    SIF input: Nick Gould, June 1994
% 
%    classification = 'OUR2-AN-V-0'
% 
%    Number of unknowns
% 
%       Alternative values for the SIF file parameters:
% IE N                   2              $-PARAMETER
% IE N                   3              $-PARAMETER
% IE N                   10             $-PARAMETER
% IE N                   100            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SENSORS';

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
% IE N                   1000           $-PARAMETER
        v_('1') = 1;
        v_('RN') = v_('N');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['THETA',int2str(I)],ix_);
            pb.xnames{iv} = ['THETA',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['S',int2str(I),',',int2str(J)],ig_);
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
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('I/N') = v_('RI')/v_('RN');
            pb.x0(ix_(['THETA',int2str(I)]),1) = v_('I/N');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSINFUN',iet_);
        elftv{it}{1} = 'THETAI';
        elftv{it}{2} = 'THETAJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                ename = ['S',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eSINFUN';
                    ielftype(ie) = iet_('eSINFUN');
                end
                vname = ['THETA',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETAI',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['THETA',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('THETAJ',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gmL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gmL2';
        end
        for J=v_('1'):v_('N')
            for I=v_('1'):v_('N')
                ig = ig_(['S',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-V-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSINFUN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TIMJ = EV_(1)-EV_(2);
        SI = sin(EV_(1));
        SJ = sin(EV_(2));
        SIMJ = sin(TIMJ);
        CI = cos(EV_(1));
        CJ = cos(EV_(2));
        CIMJ = cos(TIMJ);
        CJSIMJ = CJ*SIMJ-SJ*CIMJ;
        CJCIMJ = CJ*CIMJ+SJ*SIMJ;
        varargout{1} = SI*SJ*SIMJ;
        if(nargout>1)
            g_(1,1) = SJ*(CI*SIMJ+SI*CIMJ);
            g_(2,1) = SI*CJSIMJ;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*SJ*(CI*CIMJ-SI*SIMJ);
                H_(1,2) = CI*CJSIMJ+SI*CJCIMJ;
                H_(2,1) = H_(1,2);
                H_(2,2) = -2.0*SI*CJCIMJ;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gmL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = -GVAR_*GVAR_;
        if(nargout>1)
            g_ = -GVAR_-GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -2.0;
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

