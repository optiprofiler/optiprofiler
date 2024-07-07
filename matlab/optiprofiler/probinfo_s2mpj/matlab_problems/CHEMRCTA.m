function varargout = CHEMRCTA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHEMRCTA
%    *********
% 
%    The tubular chemical reactor model problem by Poore, using a
%    finite difference approximation to the steady state solutions.
% 
%    Source: Problem 8, eqs (8.6)--(8.9) in
%    J.J. More',
%    "A collection of nonlinear model problems"
%    Proceedings of the AMS-SIAM Summer seminar on the Computational
%    Solution of Nonlinear Systems of Equations, Colorado, 1988.
%    Argonne National Laboratory MCS-P60-0289, 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
%               minor correction by Ph. Shott, Jan 1995 and F Ruediger, Mar 1997.
% 
%    classification = 'NOR2-MN-V-V'
% 
%    The axial coordinate interval is [0,1]
% 
%    Number of discretized point for the interval [0,1].
%    The number of variables is 2N.
% 
%       Alternative values for the SIF file parameters:
% IE N                   5              $-PARAMETER n = 10
% IE N                   25             $-PARAMETER n = 50
% IE N                   50             $-PARAMETER n = 100
% IE N                   250            $-PARAMETER n = 500    original value
% IE N                   500            $-PARAMETER n = 1000
% IE N                   2500           $-PARAMETER n = 5000
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHEMRCTA';

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
        if(nargs<2)
            v_('PEM') = 1.0;  %  SIF file default value
        else
            v_('PEM') = varargin{2};
        end
        if(nargs<3)
            v_('PEH') = 5.0;  %  SIF file default value
        else
            v_('PEH') = varargin{3};
        end
        if(nargs<4)
            v_('D') = 0.135;  %  SIF file default value
        else
            v_('D') = varargin{4};
        end
        if(nargs<5)
            v_('B') = 0.5;  %  SIF file default value
        else
            v_('B') = varargin{5};
        end
        if(nargs<6)
            v_('BETA') = 2.0;  %  SIF file default value
        else
            v_('BETA') = varargin{6};
        end
        if(nargs<7)
            v_('GAMMA') = 25.0;  %  SIF file default value
        else
            v_('GAMMA') = varargin{7};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('1.0') = 1.0;
        v_('N-1') = -1+v_('N');
        v_('1/H') = v_('N-1');
        v_('-1/H') = -1.0*v_('1/H');
        v_('H') = v_('1.0')/v_('1/H');
        v_('1/H2') = v_('1/H')*v_('1/H');
        v_('-D') = -1.0*v_('D');
        v_('1/PEM') = v_('1.0')/v_('PEM');
        v_('1/H2PEM') = v_('1/PEM')*v_('1/H2');
        v_('-1/H2PM') = -1.0*v_('1/H2PEM');
        v_('HPEM') = v_('PEM')*v_('H');
        v_('-HPEM') = -1.0*v_('HPEM');
        v_('-2/H2PM') = v_('-1/H2PM')+v_('-1/H2PM');
        v_('CU1') = 1.0*v_('-HPEM');
        v_('CUI-1') = v_('1/H2PEM')+v_('1/H');
        v_('CUI') = v_('-2/H2PM')+v_('-1/H');
        v_('BD') = v_('B')*v_('D');
        v_('-BETA') = -1.0*v_('BETA');
        v_('1/PEH') = v_('1.0')/v_('PEH');
        v_('1/H2PEH') = v_('1/PEH')*v_('1/H2');
        v_('-1/H2PH') = -1.0*v_('1/H2PEH');
        v_('HPEH') = v_('PEH')*v_('H');
        v_('-HPEH') = -1.0*v_('HPEH');
        v_('-2/H2PH') = v_('-1/H2PH')+v_('-1/H2PH');
        v_('CT1') = 1.0*v_('-HPEH');
        v_('CTI-1') = v_('1/H2PEH')+v_('1/H');
        v_('CTI') = v_('-2/H2PH')+v_('-1/H');
        v_('CTI') = v_('CTI')+v_('-BETA');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['T',int2str(I)],ix_);
            pb.xnames{iv} = ['T',int2str(I)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['U',int2str(I)],ix_);
            pb.xnames{iv} = ['U',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii',['GU',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GU',int2str(round(v_('1')))];
        iv = ix_(['U',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['GU',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GU',int2str(round(v_('1')))];
        iv = ix_(['U',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('CU1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('CU1');
        end
        [ig,ig_] = s2mpjlib('ii',['GT',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GT',int2str(round(v_('1')))];
        iv = ix_(['T',int2str(round(v_('1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['GT',int2str(round(v_('1')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GT',int2str(round(v_('1')))];
        iv = ix_(['T',int2str(round(v_('2')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = v_('CT1')+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = v_('CT1');
        end
        for I=v_('2'):v_('N-1')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['GU',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['GU',int2str(I)];
            iv = ix_(['U',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CUI-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CUI-1');
            end
            iv = ix_(['U',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CUI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CUI');
            end
            iv = ix_(['U',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H2PEM')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H2PEM');
            end
            [ig,ig_] = s2mpjlib('ii',['GT',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['GT',int2str(I)];
            iv = ix_(['T',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('BETA')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('BETA');
            end
            iv = ix_(['T',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CTI-1')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CTI-1');
            end
            iv = ix_(['T',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('CTI')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('CTI');
            end
            iv = ix_(['T',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H2PEH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H2PEH');
            end
        end
        [ig,ig_] = s2mpjlib('ii',['GU',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GU',int2str(round(v_('N')))];
        iv = ix_(['U',int2str(round(v_('N-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['GU',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GU',int2str(round(v_('N')))];
        iv = ix_(['U',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['GT',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GT',int2str(round(v_('N')))];
        iv = ix_(['T',int2str(round(v_('N-1')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = -1.0;
        end
        [ig,ig_] = s2mpjlib('ii',['GT',int2str(round(v_('N')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['GT',int2str(round(v_('N')))];
        iv = ix_(['T',int2str(round(v_('N')))]);
        if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
            pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
        else
            pbm.A(ig,iv) = 1.0;
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
        pbm.gconst(ig_(['GU',int2str(round(v_('1')))])) = v_('-HPEM');
        pbm.gconst(ig_(['GT',int2str(round(v_('1')))])) = v_('-HPEH');
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['T',int2str(I)]),1) = 0.0000001;
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eREAC',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'T';
        elftp{it}{1} = 'G';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('2'):v_('N-1')
            ename = ['EU',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eREAC';
                ielftype(ie) = iet_('eREAC');
            end
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('G',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('GAMMA');
            ename = ['ET',int2str(I)];
            [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
            if(newelt)
                pbm.elftype{ie} = 'eREAC';
                ielftype(ie) = iet_('eREAC');
            end
            vname = ['U',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['T',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('G',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('GAMMA');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('2'):v_('N-1')
            ig = ig_(['GU',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EU',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-D');
            ig = ig_(['GT',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ET',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('BD');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
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

    case 'eREAC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DADT = pbm.elpar{iel_}(1)/(EV_(2)*EV_(2));
        D2ADT2 = -2.0*DADT/EV_(2);
        EX = exp(pbm.elpar{iel_}(1)-pbm.elpar{iel_}(1)/EV_(2));
        UEX = EX*EV_(1);
        varargout{1} = UEX;
        if(nargout>1)
            g_(1,1) = EX;
            g_(2,1) = UEX*DADT;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = EX*DADT;
                H_(2,1) = H_(1,2);
                H_(2,2) = UEX*(DADT*DADT+D2ADT2);
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

