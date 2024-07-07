function varargout = OSBORNEB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : OSBORNEB
%    *********
% 
%    Osborne second problem in 11 variables.
% 
%    This function  is a nonlinear least squares with 65 groups.  Each
%    group has 4 nonlinear elements.
% 
%    Source:  Problem 19 in
%    J.J. More', B.S. Garbow and K.E. Hillstrom,
%    "Testing Unconstrained Optimization Software",
%    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
%    See also Buckley#32 (p.78).
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SUR2-MN-11-0'
% 
%    Number of groups
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'OSBORNEB';

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
        v_('M') = 65;
        v_('N') = 11;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        pbm.gconst(ig_('G1')) = 1.366;
        pbm.gconst(ig_('G2')) = 1.191;
        pbm.gconst(ig_('G3')) = 1.112;
        pbm.gconst(ig_('G4')) = 1.013;
        pbm.gconst(ig_('G5')) = 0.991;
        pbm.gconst(ig_('G6')) = 0.885;
        pbm.gconst(ig_('G7')) = 0.831;
        pbm.gconst(ig_('G8')) = 0.847;
        pbm.gconst(ig_('G9')) = 0.786;
        pbm.gconst(ig_('G10')) = 0.725;
        pbm.gconst(ig_('G11')) = 0.746;
        pbm.gconst(ig_('G12')) = 0.679;
        pbm.gconst(ig_('G13')) = 0.608;
        pbm.gconst(ig_('G14')) = 0.655;
        pbm.gconst(ig_('G15')) = 0.616;
        pbm.gconst(ig_('G16')) = 0.606;
        pbm.gconst(ig_('G17')) = 0.602;
        pbm.gconst(ig_('G18')) = 0.626;
        pbm.gconst(ig_('G19')) = 0.651;
        pbm.gconst(ig_('G20')) = 0.724;
        pbm.gconst(ig_('G21')) = 0.649;
        pbm.gconst(ig_('G22')) = 0.649;
        pbm.gconst(ig_('G23')) = 0.694;
        pbm.gconst(ig_('G24')) = 0.644;
        pbm.gconst(ig_('G25')) = 0.624;
        pbm.gconst(ig_('G26')) = 0.661;
        pbm.gconst(ig_('G27')) = 0.612;
        pbm.gconst(ig_('G28')) = 0.558;
        pbm.gconst(ig_('G29')) = 0.533;
        pbm.gconst(ig_('G30')) = 0.495;
        pbm.gconst(ig_('G31')) = 0.500;
        pbm.gconst(ig_('G32')) = 0.423;
        pbm.gconst(ig_('G33')) = 0.395;
        pbm.gconst(ig_('G34')) = 0.375;
        pbm.gconst(ig_('G35')) = 0.372;
        pbm.gconst(ig_('G36')) = 0.391;
        pbm.gconst(ig_('G37')) = 0.396;
        pbm.gconst(ig_('G38')) = 0.405;
        pbm.gconst(ig_('G39')) = 0.428;
        pbm.gconst(ig_('G40')) = 0.429;
        pbm.gconst(ig_('G41')) = 0.523;
        pbm.gconst(ig_('G42')) = 0.562;
        pbm.gconst(ig_('G43')) = 0.607;
        pbm.gconst(ig_('G44')) = 0.653;
        pbm.gconst(ig_('G45')) = 0.672;
        pbm.gconst(ig_('G46')) = 0.708;
        pbm.gconst(ig_('G47')) = 0.633;
        pbm.gconst(ig_('G48')) = 0.668;
        pbm.gconst(ig_('G49')) = 0.645;
        pbm.gconst(ig_('G50')) = 0.632;
        pbm.gconst(ig_('G51')) = 0.591;
        pbm.gconst(ig_('G52')) = 0.559;
        pbm.gconst(ig_('G53')) = 0.597;
        pbm.gconst(ig_('G54')) = 0.625;
        pbm.gconst(ig_('G55')) = 0.739;
        pbm.gconst(ig_('G56')) = 0.710;
        pbm.gconst(ig_('G57')) = 0.729;
        pbm.gconst(ig_('G58')) = 0.720;
        pbm.gconst(ig_('G59')) = 0.636;
        pbm.gconst(ig_('G60')) = 0.581;
        pbm.gconst(ig_('G61')) = 0.428;
        pbm.gconst(ig_('G62')) = 0.292;
        pbm.gconst(ig_('G63')) = 0.162;
        pbm.gconst(ig_('G64')) = 0.098;
        pbm.gconst(ig_('G65')) = 0.054;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 1.3;
        pb.x0(ix_('X2'),1) = 0.65;
        pb.x0(ix_('X3'),1) = 0.65;
        pb.x0(ix_('X4'),1) = 0.7;
        pb.x0(ix_('X5'),1) = 0.6;
        pb.x0(ix_('X6'),1) = 3.0;
        pb.x0(ix_('X7'),1) = 5.0;
        pb.x0(ix_('X8'),1) = 7.0;
        pb.x0(ix_('X9'),1) = 2.0;
        pb.x0(ix_('X10'),1) = 4.5;
        pb.x0(ix_('X11'),1) = 5.5;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePEXP',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'ePEXP3',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftp{it}{1} = 'T3';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            v_('I-1') = 1+I;
            v_('RI-1') = v_('I-1');
            v_('TI') = 0.1*v_('RI-1');
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP';
            ielftype(ie) = iet_('ePEXP');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP3';
            ielftype(ie) = iet_('ePEXP3');
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X9';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP3';
            ielftype(ie) = iet_('ePEXP3');
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X10';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
            ename = ['D',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePEXP3';
            ielftype(ie) = iet_('ePEXP3');
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X11';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X8';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('TI');
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
        for I=v_('1'):v_('M')
            ig = ig_(['G',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['D',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.04013774
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-11-0';
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

    case 'ePEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPA = exp(-pbm.elpar{iel_}(1)*EV_(2));
        FVAL = EV_(1)*EXPA;
        varargout{1} = FVAL;
        if(nargout>1)
            g_(1,1) = EXPA;
            g_(2,1) = -pbm.elpar{iel_}(1)*FVAL;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -pbm.elpar{iel_}(1)*EXPA;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*FVAL;
                varargout{3} = H_;
            end
        end

    case 'ePEXP3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TMV2 = pbm.elpar{iel_}(1)-EV_(2);
        TMV2SQ = TMV2*TMV2;
        EXPA = exp(-TMV2SQ*EV_(3));
        FVAL = EV_(1)*EXPA;
        A = 2.0*TMV2*EV_(3);
        varargout{1} = FVAL;
        if(nargout>1)
            g_(1,1) = EXPA;
            g_(2,1) = A*FVAL;
            g_(3,1) = -TMV2SQ*FVAL;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = A*EXPA;
                H_(2,1) = H_(1,2);
                H_(1,3) = -TMV2SQ*EXPA;
                H_(3,1) = H_(1,3);
                H_(2,2) = (A*A-2.0*EV_(3))*FVAL;
                H_(2,3) = (2.0*TMV2-A*TMV2SQ)*FVAL;
                H_(3,2) = H_(2,3);
                H_(3,3) = TMV2SQ*TMV2SQ*FVAL;
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

