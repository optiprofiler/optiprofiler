function varargout = HIMMELBI(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELBI
%    *********
% 
%    An unpleasant weapon assignment problem by Bracken and McCormick.
% 
%    The real problem has integer variables.
%    Also, the sign of ci have been reversed in order to have a
%    meaningful constraints on the total number of weapons (a fully
%    desirable situation).
% 
%    Source: problem 23 in
%    D.H. Himmelblau,
%    "Applied nonlinear programming",
%    McGraw-Hill, New-York, 1972.
% 
%    SIF input: Ph. Toint, March 1991.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-COLR2-MN-100-12'
% 
%    Number of variables
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELBI';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('N') = 100;
        v_('NT') = 20;
        v_('A1,1') = 1.0;
        v_('A2,1') = 0.84;
        v_('A3,1') = 0.96;
        v_('A4,1') = 1.0;
        v_('A5,1') = 0.92;
        v_('A1,2') = 0.95;
        v_('A2,2') = 0.83;
        v_('A3,2') = 0.95;
        v_('A4,2') = 1.0;
        v_('A5,2') = 0.94;
        v_('A1,3') = 1.0;
        v_('A2,3') = 0.85;
        v_('A3,3') = 0.95;
        v_('A4,3') = 1.0;
        v_('A5,3') = 0.92;
        v_('A1,4') = 1.0;
        v_('A2,4') = 0.84;
        v_('A3,4') = 0.96;
        v_('A4,4') = 1.0;
        v_('A5,4') = 0.95;
        v_('A1,5') = 1.0;
        v_('A2,5') = 0.85;
        v_('A3,5') = 0.96;
        v_('A4,5') = 1.0;
        v_('A5,5') = 0.95;
        v_('A1,6') = 0.85;
        v_('A2,6') = 0.81;
        v_('A3,6') = 0.90;
        v_('A4,6') = 1.0;
        v_('A5,6') = 0.98;
        v_('A1,7') = 0.90;
        v_('A2,7') = 0.81;
        v_('A3,7') = 0.92;
        v_('A4,7') = 1.0;
        v_('A5,7') = 0.98;
        v_('A1,8') = 0.85;
        v_('A2,8') = 0.82;
        v_('A3,8') = 0.91;
        v_('A4,8') = 1.0;
        v_('A5,8') = 1.0;
        v_('A1,9') = 0.80;
        v_('A2,9') = 0.80;
        v_('A3,9') = 0.92;
        v_('A4,9') = 1.0;
        v_('A5,9') = 1.0;
        v_('A1,10') = 1.0;
        v_('A2,10') = 0.86;
        v_('A3,10') = 0.95;
        v_('A4,10') = 0.96;
        v_('A5,10') = 0.90;
        v_('A1,11') = 1.0;
        v_('A2,11') = 1.0;
        v_('A3,11') = 0.99;
        v_('A4,11') = 0.91;
        v_('A5,11') = 0.95;
        v_('A1,12') = 1.0;
        v_('A2,12') = 0.98;
        v_('A3,12') = 0.98;
        v_('A4,12') = 0.92;
        v_('A5,12') = 0.96;
        v_('A1,13') = 1.0;
        v_('A2,13') = 1.0;
        v_('A3,13') = 0.99;
        v_('A4,13') = 0.91;
        v_('A5,13') = 0.91;
        v_('A1,14') = 1.0;
        v_('A2,14') = 0.88;
        v_('A3,14') = 0.98;
        v_('A4,14') = 0.92;
        v_('A5,14') = 0.98;
        v_('A1,15') = 1.0;
        v_('A2,15') = 0.87;
        v_('A3,15') = 0.97;
        v_('A4,15') = 0.98;
        v_('A5,15') = 0.99;
        v_('A1,16') = 1.0;
        v_('A2,16') = 0.88;
        v_('A3,16') = 0.98;
        v_('A4,16') = 0.93;
        v_('A5,16') = 0.99;
        v_('A1,17') = 1.0;
        v_('A2,17') = 0.85;
        v_('A3,17') = 0.95;
        v_('A4,17') = 1.0;
        v_('A5,17') = 1.0;
        v_('A1,18') = 0.95;
        v_('A2,18') = 0.84;
        v_('A3,18') = 0.92;
        v_('A4,18') = 1.0;
        v_('A5,18') = 1.0;
        v_('A1,19') = 1.0;
        v_('A2,19') = 0.85;
        v_('A3,19') = 0.93;
        v_('A4,19') = 1.0;
        v_('A5,19') = 1.0;
        v_('A1,20') = 1.0;
        v_('A2,20') = 0.85;
        v_('A3,20') = 0.92;
        v_('A4,20') = 1.0;
        v_('A5,20') = 1.0;
        v_('B1') = 30.0;
        v_('B6') = 100.0;
        v_('B10') = 40.0;
        v_('B14') = 50.0;
        v_('B15') = 70.0;
        v_('B16') = 35.0;
        v_('B20') = 10.0;
        v_('U1') = 60.0;
        v_('U2') = 50.0;
        v_('U3') = 50.0;
        v_('U4') = 75.0;
        v_('U5') = 40.0;
        v_('U6') = 60.0;
        v_('U7') = 35.0;
        v_('U8') = 30.0;
        v_('U9') = 25.0;
        v_('U10') = 150.0;
        v_('U11') = 30.0;
        v_('U12') = 45.0;
        v_('U13') = 125.0;
        v_('U14') = 200.0;
        v_('U15') = 200.0;
        v_('U16') = 130.0;
        v_('U17') = 100.0;
        v_('U18') = 100.0;
        v_('U19') = 100.0;
        v_('U20') = 150.0;
        v_('C1') = 200.0;
        v_('C2') = 100.0;
        v_('C3') = 300.0;
        v_('C4') = 150.0;
        v_('C5') = 250.0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('NW') = 0.0;
        for I=v_('1'):v_('5')
            v_('NW') = v_('NW')+v_(['C',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for J=v_('1'):v_('NT')
            for I=v_('1'):v_('5')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for J=v_('1'):v_('NT')
            [ig,ig_] = s2mpjlib('ii',['P',int2str(J)],ig_);
            gtype{ig} = '<>';
            v_('1/UJ') = 1.0/v_(['U',int2str(J)]);
            pbm.gscale(ig,1) = v_('1/UJ');
        end
        v_('J') = 1;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 6;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 10;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 14;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 15;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 16;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        v_('J') = 20;
        for I=v_('1'):v_('5')
            [ig,ig_] = s2mpjlib('ii',['CB',int2str(round(v_('J')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['CB',int2str(round(v_('J')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I),',',int2str(round(v_('J')))]);
            valA(end+1) = 1.0;
        end
        for I=v_('1'):v_('5')
            for J=v_('1'):v_('NT')
                [ig,ig_] = s2mpjlib('ii',['CC',int2str(I)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['CC',int2str(I)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['X',int2str(I),',',int2str(J)]);
                valA(end+1) = 1.0;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        for J=v_('1'):v_('NT')
            pbm.gconst(ig_(['P',int2str(J)])) = 1.0;
        end
        v_('J') = 1;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 6;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 10;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 14;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 15;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 16;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        v_('J') = 20;
        pbm.gconst(ig_(['CB',int2str(round(v_('J')))])) =...
              v_(['B',int2str(round(v_('J')))]);
        for I=v_('1'):v_('5')
            pbm.gconst(ig_(['CC',int2str(I)])) = v_(['C',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = v_('NW')*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en5PEXP',iet_);
        elftv{it}{1} = 'Y1';
        elftv{it}{2} = 'Y2';
        elftv{it}{3} = 'Y3';
        elftv{it}{4} = 'Y4';
        elftv{it}{5} = 'Y5';
        elftp{it}{1} = 'A1';
        elftp{it}{2} = 'A2';
        elftp{it}{3} = 'A3';
        elftp{it}{4} = 'A4';
        elftp{it}{5} = 'A5';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for J=v_('1'):v_('NT')
            ename = ['PP',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en5PEXP';
            ielftype(ie) = iet_('en5PEXP');
            vname = ['X',int2str(round(v_('1'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],v_('NW'),[]);
            posev = find(strcmp('Y1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('2'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],v_('NW'),[]);
            posev = find(strcmp('Y2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('3'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],v_('NW'),[]);
            posev = find(strcmp('Y3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('4'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],v_('NW'),[]);
            posev = find(strcmp('Y4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('5'))),',',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],v_('NW'),[]);
            posev = find(strcmp('Y5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(round(v_('1'))),',',int2str(J)]);
            [~,posep] = ismember('A2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(round(v_('2'))),',',int2str(J)]);
            [~,posep] = ismember('A3',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(round(v_('3'))),',',int2str(J)]);
            [~,posep] = ismember('A4',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(round(v_('4'))),',',int2str(J)]);
            [~,posep] = ismember('A5',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['A',int2str(round(v_('5'))),',',int2str(J)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('1'):v_('NT')
            ig = ig_(['P',int2str(J)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['PP',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               -1735.56958
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-MN-100-12';
        pb.x0          = zeros(pb.n,1);
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en5PEXP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LA1 = log(pbm.elpar{iel_}(1));
        LA2 = log(pbm.elpar{iel_}(2));
        LA3 = log(pbm.elpar{iel_}(3));
        LA4 = log(pbm.elpar{iel_}(4));
        LA5 = log(pbm.elpar{iel_}(5));
        FF =...
              pbm.elpar{iel_}(1)^EV_(1)*pbm.elpar{iel_}(2)^EV_(2)*pbm.elpar{iel_}(3)^EV_(3)*pbm.elpar{iel_}(4)^EV_(4)*pbm.elpar{iel_}(5)^EV_(5);
        varargout{1} = FF;
        if(nargout>1)
            g_(1,1) = LA1*FF;
            g_(2,1) = LA2*FF;
            g_(3,1) = LA3*FF;
            g_(4,1) = LA4*FF;
            g_(5,1) = LA5*FF;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = LA1*LA1*FF;
                H_(1,2) = LA1*LA2*FF;
                H_(2,1) = H_(1,2);
                H_(1,3) = LA1*LA3*FF;
                H_(3,1) = H_(1,3);
                H_(1,4) = LA1*LA4*FF;
                H_(4,1) = H_(1,4);
                H_(1,5) = LA1*LA5*FF;
                H_(5,1) = H_(1,5);
                H_(2,2) = LA2*LA2*FF;
                H_(2,3) = LA2*LA3*FF;
                H_(3,2) = H_(2,3);
                H_(2,4) = LA2*LA4*FF;
                H_(4,2) = H_(2,4);
                H_(2,5) = LA2*LA5*FF;
                H_(5,2) = H_(2,5);
                H_(3,3) = LA3*LA3*FF;
                H_(3,4) = LA3*LA4*FF;
                H_(4,3) = H_(3,4);
                H_(3,5) = LA3*LA5*FF;
                H_(5,3) = H_(3,5);
                H_(4,4) = LA4*LA4*FF;
                H_(4,5) = LA4*LA5*FF;
                H_(5,4) = H_(4,5);
                H_(5,5) = LA5*LA5*FF;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

