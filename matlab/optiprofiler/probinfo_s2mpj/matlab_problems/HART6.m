function varargout = HART6(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HART6
%    *********
% 
%    Source: Hartman problem 6 in
%    L. C. W. Dixon and G. P. Szego (Eds.)
%    Towards Global Optimization
%    North Holland, 1975.
%    Paper 9, page 163.
% 
%    SIF input: A.R. Conn May 1995
% 
%    classification = 'OBR2-AN-6-0'
% 
%    Number of variables - constraints
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HART6';

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
        v_('N') = 6;
        v_('1') = 1;
        v_('ONE') = 1;
        v_('NN') = 6;
        v_('L') = 4;
        v_('C1') = 1.0;
        v_('C2') = 1.2;
        v_('C3') = 3.0;
        v_('C4') = 3.2;
        v_('A1,1') = 10.0;
        v_('A2,1') = 0.05;
        v_('A3,1') = 3.0;
        v_('A4,1') = 17.0;
        v_('A1,2') = 0.05;
        v_('A2,2') = 10.0;
        v_('A3,2') = 3.5;
        v_('A4,2') = 8.0;
        v_('A1,3') = 17.0;
        v_('A2,3') = 17.0;
        v_('A3,3') = 1.7;
        v_('A4,3') = 0.05;
        v_('A1,4') = 3.5;
        v_('A2,4') = 0.1;
        v_('A3,4') = 10.0;
        v_('A4,4') = 10.0;
        v_('A1,5') = 1.7;
        v_('A2,5') = 8.0;
        v_('A3,5') = 17.0;
        v_('A4,5') = 0.1;
        v_('A1,6') = 8.0;
        v_('A2,6') = 14.0;
        v_('A3,6') = 8.0;
        v_('A4,6') = 14.0;
        v_('P1,1') = 0.1312;
        v_('P2,1') = 0.2329;
        v_('P3,1') = 0.2348;
        v_('P4,1') = 0.4047;
        v_('P1,2') = 0.1696;
        v_('P2,2') = 0.4135;
        v_('P3,2') = 0.1451;
        v_('P4,2') = 0.8828;
        v_('P1,3') = 0.5569;
        v_('P2,3') = 0.8307;
        v_('P3,3') = 0.3522;
        v_('P4,3') = 0.8732;
        v_('P1,4') = 0.0124;
        v_('P2,4') = 0.3736;
        v_('P3,4') = 0.2883;
        v_('P4,4') = 0.5743;
        v_('P1,5') = 0.8283;
        v_('P2,5') = 0.1004;
        v_('P3,5') = 0.3047;
        v_('P4,5') = 0.1091;
        v_('P1,6') = 0.5886;
        v_('P2,6') = 0.9991;
        v_('P3,6') = 0.6650;
        v_('P4,6') = 0.0381;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('L')
            [ig,ig_] = s2mpjlib('ii',['OBJ',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = -1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.2*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'V1';
        elftp{it}{1} = 'PIJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('L')
            for J=v_('1'):v_('NN')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,0.0,1.0,0.2);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('PIJ',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['P',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gNEXP',igt_);
        [it,igt_] = s2mpjlib('ii','gNEXP',igt_);
        grftp{it}{1} = 'CI';
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('L')
            ig = ig_(['OBJ',int2str(I)]);
            pbm.grftype{ig} = 'gNEXP';
            for J=v_('1'):v_('NN')
                ig = ig_(['OBJ',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = v_(['A',int2str(I),',',int2str(J)]);
            end
            ig = ig_(['OBJ',int2str(I)]);
            [~,posgp] = ismember('CI',grftp{igt_(pbm.grftype{ig})});
            pbm.grpar{ig}(posgp) = v_(['C',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -3.32288689158
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-AN-6-0';
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
        varargout{1} = (EV_(1)-pbm.elpar{iel_}(1))*(EV_(1)-pbm.elpar{iel_}(1));
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-pbm.elpar{iel_}(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gNEXP'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = pbm.grpar{igr_}(1)*exp(-GVAR_);
        if(nargout>1)
            g_ = -pbm.grpar{igr_}(1)*exp(-GVAR_);
            varargout{2} = g_;
            if(nargout>2)
                H_ = pbm.grpar{igr_}(1)*exp(-GVAR_);
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

