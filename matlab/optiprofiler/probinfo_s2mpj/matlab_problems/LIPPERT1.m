function varargout = LIPPERT1(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LIPPERT1
%    *********
% 
%    A discrete approximation to a continuum optimal flow problem
%    in the unit square. The continuum problem requires that the
%    divergence of a given flow should be given everywhere in the
%    region of interest, with the restriction that the capacity of
%    the flow is bounded. The aim is then to maximize the given flow.
% 
%    The discrete problem (primal formulation 1) in the unit square is to 
%      maximize   t
%      subject to dx( u_ij - ui-1j ) + dx( v_ij - vij-1 ) = t s_ij
%                 u_ij^2 + v_ij^2 <= 1
%                 u_i-1j^2 + v_ij^2 <= 1
%                 u_ij^2 + v_ij-1^2 <= 1
%                 u_i-1j^2 + v_ij-1^2 <= 1
%      where 1 <= i <= nx, 1 <= j <= ny
%      and        t >= 0
% 
%    Source: R. A. Lippert
%      "Discrete approximations to continuum optimal flow problems"
%      Tech. Report, Dept of Maths, M.I.T., 2006
%    following a suggestion by Gil Strang
% 
%    SIF input: Nick Gould, September 2006
% 
%    classification = 'LQR2-MN-V-V'
% 
%    Number of nodes in x direction
% 
%       Alternative values for the SIF file parameters:
% IE NX                  2              $-PARAMETER
% IE NX                  3              $-PARAMETER
% IE NX                  10             $-PARAMETER
% IE NX                  40             $-PARAMETER
% IE NX                  100            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LIPPERT1';

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
            v_('NX') = 3;  %  SIF file default value
        else
            v_('NX') = varargin{1};
        end
% IE NY                  2              $-PARAMETER
% IE NY                  3              $-PARAMETER
% IE NY                  10             $-PARAMETER 
% IE NY                  40             $-PARAMETER
% IE NY                  100            $-PARAMETER
        if(nargs<2)
            v_('NY') = 10;  %  SIF file default value
        else
            v_('NY') = varargin{2};
        end
        v_('X+') = 1+v_('NX');
        v_('X-') = -1+v_('NX');
        v_('Y+') = 1+v_('NY');
        v_('Y-') = -1+v_('NY');
        v_('1') = 1;
        v_('0') = 0;
        v_('ONE') = 1.0;
        v_('-ONE') = -1.0;
        v_('S') = 1.0;
        v_('-S') = v_('S')*v_('-ONE');
        v_('RX') = v_('NX');
        v_('DX') = v_('ONE')/v_('RX');
        v_('-DX') = v_('-ONE')/v_('RX');
        v_('RY') = v_('NY');
        v_('DY') = v_('ONE')/v_('RY');
        v_('-DY') = v_('-ONE')/v_('RY');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','T',ix_);
        pb.xnames{iv} = 'T';
        for I=v_('0'):v_('NX')
            for J=v_('1'):v_('NY')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('NX')
            for J=v_('0'):v_('NY')
                [iv,ix_] = s2mpjlib('ii',['V',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['V',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        iv = ix_('T');
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        for I=v_('1'):v_('NX')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('NY')
                v_('J-1') = -1+J;
                [ig,ig_] = s2mpjlib('ii',['O',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['O',int2str(I),',',int2str(J)];
                pbm.gscale(ig,1) = v_('DX');
                iv = ix_(['U',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('DX');
                end
                iv = ix_(['U',int2str(round(v_('I-1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-DX');
                end
                iv = ix_(['V',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('DY');
                end
                iv = ix_(['V',int2str(I),',',int2str(round(v_('J-1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-DY');
                end
                iv = ix_('T');
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-S')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-S');
                end
            end
        end
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['A',int2str(I),',',int2str(J)];
                [ig,ig_] = s2mpjlib('ii',['B',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['B',int2str(I),',',int2str(J)];
                [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['C',int2str(I),',',int2str(J)];
                [ig,ig_] = s2mpjlib('ii',['D',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['D',int2str(I),',',int2str(J)];
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
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('NX')
            for J=v_('1'):v_('NY')
                pbm.gconst(ig_(['A',int2str(I),',',int2str(J)])) = 1.0;
                pbm.gconst(ig_(['B',int2str(I),',',int2str(J)])) = 1.0;
                pbm.gconst(ig_(['C',int2str(I),',',int2str(J)])) = 1.0;
                pbm.gconst(ig_(['D',int2str(I),',',int2str(J)])) = 1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_('T'),1) = 0.01;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'ALPHA';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('0'):v_('NX')
            for J=v_('1'):v_('NY')
                ename = ['P',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eSQR';
                    ielftype(ie) = iet_('eSQR');
                end
                vname = ['U',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ALPHA',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('1'):v_('NX')
            for J=v_('0'):v_('NY')
                ename = ['Q',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'eSQR';
                    ielftype(ie) = iet_('eSQR');
                end
                vname = ['V',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('ALPHA',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NX')
            v_('I-1') = -1+I;
            for J=v_('1'):v_('NY')
                v_('J-1') = -1+J;
                ig = ig_(['A',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['B',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I-1'))),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['C',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(round(v_('J-1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                ig = ig_(['D',int2str(I),',',int2str(J)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I-1'))),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['Q',int2str(I),',',int2str(round(v_('J-1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -3.77245385
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
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

