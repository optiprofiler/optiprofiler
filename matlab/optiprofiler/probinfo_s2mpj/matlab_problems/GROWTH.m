function varargout = GROWTH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : GROWTH
%    *********
%    GROWTH problem in 3 variables
% 
%    Fit the observed growth g(n) from Gaussian Elimination
%    with complete pivoting to a function of the form
%         U1 * n ** ( U2 + LOG(n) * U3 )
% 
%    SIF input: Nick Gould, Nov, 1991.
% 
%    classification = 'C-CNOR2-AN-3-12'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'GROWTH';

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
        v_('N') = 3;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','U1',ix_);
        pb.xnames{iv} = 'U1';
        [iv,ix_] = s2mpjlib('ii','U2',ix_);
        pb.xnames{iv} = 'U2';
        [iv,ix_] = s2mpjlib('ii','U3',ix_);
        pb.xnames{iv} = 'U3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','G8',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G8';
        [ig,ig_] = s2mpjlib('ii','G9',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G9';
        [ig,ig_] = s2mpjlib('ii','G10',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G10';
        [ig,ig_] = s2mpjlib('ii','G11',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G11';
        [ig,ig_] = s2mpjlib('ii','G12',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G12';
        [ig,ig_] = s2mpjlib('ii','G13',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G13';
        [ig,ig_] = s2mpjlib('ii','G14',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G14';
        [ig,ig_] = s2mpjlib('ii','G15',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G15';
        [ig,ig_] = s2mpjlib('ii','G16',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G16';
        [ig,ig_] = s2mpjlib('ii','G18',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G18';
        [ig,ig_] = s2mpjlib('ii','G20',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G20';
        [ig,ig_] = s2mpjlib('ii','G25',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'G25';
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
        pbm.gconst(ig_('G8')) = 8.0;
        pbm.gconst(ig_('G9')) = 8.4305;
        pbm.gconst(ig_('G10')) = 9.5294;
        pbm.gconst(ig_('G11')) = 10.4627;
        pbm.gconst(ig_('G12')) = 12.0;
        pbm.gconst(ig_('G13')) = 13.0205;
        pbm.gconst(ig_('G14')) = 14.5949;
        pbm.gconst(ig_('G15')) = 16.1078;
        pbm.gconst(ig_('G16')) = 18.0596;
        pbm.gconst(ig_('G18')) = 20.4569;
        pbm.gconst(ig_('G20')) = 24.25;
        pbm.gconst(ig_('G25')) = 32.9863;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('U1'),1) = 100.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eFIT',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftv{it}{3} = 'U3';
        elftp{it}{1} = 'RN';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'G8';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 8.0;
        ename = 'G9';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 9.0;
        ename = 'G10';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 10.0;
        ename = 'G11';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 11.0;
        ename = 'G12';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 12.0;
        ename = 'G13';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 13.0;
        ename = 'G14';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 14.0;
        ename = 'G15';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 15.0;
        ename = 'G16';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 16.0;
        ename = 'G18';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 18.0;
        ename = 'G20';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 20.0;
        ename = 'G25';
        [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
        if(newelt)
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
        end
        vname = 'U1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U1',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U2',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'U3';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('U3',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('RN',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = 25.0;
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        ig = ig_('G8');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G8');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G9');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G9');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G10');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G10');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G11');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G11');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G12');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G12');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G13');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G13');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G14');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G14');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G15');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G15');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G16');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G16');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G18');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G18');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G20');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G20');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_('G25');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('G25');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-AN-3-12';
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

    case 'eFIT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        LOGRN = log(pbm.elpar{iel_}(1));
        POWER = pbm.elpar{iel_}(1)^(EV_(2)+LOGRN*EV_(3));
        varargout{1} = EV_(1)*POWER;
        if(nargout>1)
            g_(1,1) = POWER;
            g_(2,1) = EV_(1)*POWER*LOGRN;
            g_(3,1) = EV_(1)*POWER*LOGRN^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = 0.0;
                H_(1,2) = POWER*LOGRN;
                H_(2,1) = H_(1,2);
                H_(1,3) = POWER*LOGRN^2;
                H_(3,1) = H_(1,3);
                H_(2,2) = EV_(1)*POWER*LOGRN^2;
                H_(2,3) = EV_(1)*POWER*LOGRN^3;
                H_(3,2) = H_(2,3);
                H_(3,3) = EV_(1)*POWER*LOGRN^4;
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

