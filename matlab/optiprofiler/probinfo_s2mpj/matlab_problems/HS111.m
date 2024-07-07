function varargout = HS111(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS111
%    *********
% 
%    This problem is a chemical equilibrium problem involving 3 linear
%    equality constraints.
% 
%    Source: problem 111 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    SIF input: Nick Gould, August 1991.
% 
%    classification = 'OOR2-AN-10-3'
% 
%    N is the number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS111';

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
        v_('N') = 10;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('C1') = -6.089;
        v_('C2') = -17.164;
        v_('C3') = -34.054;
        v_('C4') = -5.914;
        v_('C5') = -24.721;
        v_('C6') = -14.986;
        v_('C7') = -24.100;
        v_('C8') = -10.708;
        v_('C9') = -26.662;
        v_('C10') = -22.179;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        [ig,ig_] = s2mpjlib('ii','CON1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON1';
        [ig,ig_] = s2mpjlib('ii','CON2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON2';
        [ig,ig_] = s2mpjlib('ii','CON3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'CON3';
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
        pbm.gconst(ig_('CON1')) = 2.0;
        pbm.gconst(ig_('CON2')) = 1.0;
        pbm.gconst(ig_('CON3')) = 1.0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -100.0*ones(pb.n,1);
        pb.xupper = 100.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = -2.3*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eOBJ',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftv{it}{5} = 'V5';
        elftv{it}{6} = 'V6';
        elftv{it}{7} = 'V7';
        elftv{it}{8} = 'V8';
        elftv{it}{9} = 'V9';
        elftv{it}{10} = 'V10';
        elftp{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eEXP',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('N')
            v_(['RI',int2str(I)]) = I;
            v_(['RI',int2str(I)]) = 0.1+v_(['RI',int2str(I)]);
        end
        for I=v_('1'):v_('N')
            ename = ['O',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eOBJ';
            ielftype(ie) = iet_('eOBJ');
            v_('TEMP') = v_(['RI',int2str(round(v_('1')))]);
            v_(['RI',int2str(round(v_('1')))]) = v_(['RI',int2str(I)]);
            v_(['RI',int2str(I)]) = v_('TEMP');
            v_('R') = v_(['RI',int2str(round(v_('1')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('2')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('3')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('4')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('5')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('6')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('7')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V7',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('8')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V8',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('9')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V9',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('R') = v_(['RI',int2str(round(v_('10')))]);
            v_('J') = fix(v_('R'));
            vname = ['X',int2str(round(v_('J')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('V10',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('C',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['C',int2str(I)]);
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXP';
            ielftype(ie) = iet_('eEXP');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-100.0,100.0,-2.3);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('N')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['O',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_('CON1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E2');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 2.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('CON2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E5');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        pbm.grelw{ig}(posel) = 1.0;
        ig = ig_('CON3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E7');
        pbm.grelw{ig}(posel) = 1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('E9');
        pbm.grelw{ig}(posel) = 2.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -47.707579
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-AN-10-3';
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

    case 'eEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EX = exp(EV_(1));
        varargout{1} = EX;
        if(nargout>1)
            g_(1,1) = EX;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = EX;
                varargout{3} = H_;
            end
        end

    case 'eOBJ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E1 = exp(EV_(1));
        E2 = exp(EV_(2));
        E3 = exp(EV_(3));
        E4 = exp(EV_(4));
        E5 = exp(EV_(5));
        E6 = exp(EV_(6));
        E7 = exp(EV_(7));
        E8 = exp(EV_(8));
        E9 = exp(EV_(9));
        E10 = exp(EV_(10));
        SUM = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10;
        varargout{1} = E1*(pbm.elpar{iel_}(1)+EV_(1)-log(SUM));
        if(nargout>1)
            g_(1,1) = E1*(pbm.elpar{iel_}(1)+EV_(1)-log(SUM))+E1*(1.0e+0-E1/SUM);
            g_(2,1) = -E1*E2/SUM;
            g_(3,1) = -E1*E3/SUM;
            g_(4,1) = -E1*E4/SUM;
            g_(5,1) = -E1*E5/SUM;
            g_(6,1) = -E1*E6/SUM;
            g_(7,1) = -E1*E7/SUM;
            g_(8,1) = -E1*E8/SUM;
            g_(9,1) = -E1*E9/SUM;
            g_(10,1) = -E1*E10/SUM;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(10,10);
                H_(1,1) = E1*(pbm.elpar{iel_}(1)+EV_(1)-log(SUM))+E1*(1.0e+0-E1/SUM)+...
                     E1*(1.0e+0-E1/SUM)+E1*(-E1/SUM)+E1*(E1^2/SUM^2);
                H_(1,2) = (-1.0e+0+E1/SUM)*E1*E2/SUM;
                H_(2,1) = H_(1,2);
                H_(2,2) = (-1.0e+0+E2/SUM)*E1*E2/SUM;
                H_(1,3) = (-1.0e+0+E1/SUM)*E1*E3/SUM;
                H_(3,1) = H_(1,3);
                H_(2,3) = E1*E2*E3/SUM^2;
                H_(3,2) = H_(2,3);
                H_(3,3) = (-1.0e+0+E3/SUM)*E1*E3/SUM;
                H_(1,4) = (-1.0e+0+E1/SUM)*E1*E4/SUM;
                H_(4,1) = H_(1,4);
                H_(2,4) = E1*E2*E4/SUM^2;
                H_(4,2) = H_(2,4);
                H_(3,4) = E1*E3*E4/SUM^2;
                H_(4,3) = H_(3,4);
                H_(4,4) = (-1.0e+0+E4/SUM)*E1*E4/SUM;
                H_(1,5) = (-1.0e+0+E1/SUM)*E1*E5/SUM;
                H_(5,1) = H_(1,5);
                H_(2,5) = E1*E2*E5/SUM^2;
                H_(5,2) = H_(2,5);
                H_(3,5) = E1*E3*E5/SUM^2;
                H_(5,3) = H_(3,5);
                H_(4,5) = E1*E4*E5/SUM^2;
                H_(5,4) = H_(4,5);
                H_(5,5) = (-1.0e+0+E5/SUM)*E1*E5/SUM;
                H_(1,6) = (-1.0e+0+E1/SUM)*E1*E6/SUM;
                H_(6,1) = H_(1,6);
                H_(2,6) = E1*E2*E6/SUM^2;
                H_(6,2) = H_(2,6);
                H_(3,6) = E1*E3*E6/SUM^2;
                H_(6,3) = H_(3,6);
                H_(4,6) = E1*E4*E6/SUM^2;
                H_(6,4) = H_(4,6);
                H_(5,6) = E1*E5*E6/SUM^2;
                H_(6,5) = H_(5,6);
                H_(6,6) = (-1.0e+0+E6/SUM)*E1*E6/SUM;
                H_(1,7) = (-1.0e+0+E1/SUM)*E1*E7/SUM;
                H_(7,1) = H_(1,7);
                H_(2,7) = E1*E2*E7/SUM^2;
                H_(7,2) = H_(2,7);
                H_(3,7) = E1*E3*E7/SUM^2;
                H_(7,3) = H_(3,7);
                H_(4,7) = E1*E4*E7/SUM^2;
                H_(7,4) = H_(4,7);
                H_(5,7) = E1*E5*E7/SUM^2;
                H_(7,5) = H_(5,7);
                H_(6,7) = E1*E6*E7/SUM^2;
                H_(7,6) = H_(6,7);
                H_(7,7) = (-1.0e+0+E7/SUM)*E1*E7/SUM;
                H_(1,8) = (-1.0e+0+E1/SUM)*E1*E8/SUM;
                H_(8,1) = H_(1,8);
                H_(2,8) = E1*E2*E8/SUM^2;
                H_(8,2) = H_(2,8);
                H_(3,8) = E1*E3*E8/SUM^2;
                H_(8,3) = H_(3,8);
                H_(4,8) = E1*E4*E8/SUM^2;
                H_(8,4) = H_(4,8);
                H_(5,8) = E1*E5*E8/SUM^2;
                H_(8,5) = H_(5,8);
                H_(6,8) = E1*E6*E8/SUM^2;
                H_(8,6) = H_(6,8);
                H_(7,8) = E1*E7*E8/SUM^2;
                H_(8,7) = H_(7,8);
                H_(8,8) = (-1.0e+0+E8/SUM)*E1*E8/SUM;
                H_(1,9) = (-1.0e+0+E1/SUM)*E1*E9/SUM;
                H_(9,1) = H_(1,9);
                H_(2,9) = E1*E2*E9/SUM^2;
                H_(9,2) = H_(2,9);
                H_(3,9) = E1*E3*E9/SUM^2;
                H_(9,3) = H_(3,9);
                H_(4,9) = E1*E4*E9/SUM^2;
                H_(9,4) = H_(4,9);
                H_(5,9) = E1*E5*E9/SUM^2;
                H_(9,5) = H_(5,9);
                H_(6,9) = E1*E6*E9/SUM^2;
                H_(9,6) = H_(6,9);
                H_(7,9) = E1*E7*E9/SUM^2;
                H_(9,7) = H_(7,9);
                H_(8,9) = E1*E8*E9/SUM^2;
                H_(9,8) = H_(8,9);
                H_(9,9) = (-1.0e+0+E9/SUM)*E1*E9/SUM;
                H_(1,10) = (-1.0e+0+E1/SUM)*E1*E10/SUM;
                H_(10,1) = H_(1,10);
                H_(2,10) = E1*E2*E10/SUM^2;
                H_(10,2) = H_(2,10);
                H_(3,10) = E1*E3*E10/SUM^2;
                H_(10,3) = H_(3,10);
                H_(4,10) = E1*E4*E10/SUM^2;
                H_(10,4) = H_(4,10);
                H_(5,10) = E1*E5*E10/SUM^2;
                H_(10,5) = H_(5,10);
                H_(6,10) = E1*E6*E10/SUM^2;
                H_(10,6) = H_(6,10);
                H_(7,10) = E1*E7*E10/SUM^2;
                H_(10,7) = H_(7,10);
                H_(8,10) = E1*E8*E10/SUM^2;
                H_(10,8) = H_(8,10);
                H_(9,10) = E1*E9*E10/SUM^2;
                H_(10,9) = H_(9,10);
                H_(10,10) = (-1.0e+0+E10/SUM)*E1*E10/SUM;
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

