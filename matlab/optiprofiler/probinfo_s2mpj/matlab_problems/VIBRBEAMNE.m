function varargout = VIBRBEAMNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem:
%    ********
% 
%    A nonlinear least-squares problem arising from laser-Doppler
%    measurements of a vibrating beam.  The data correspond to a simulated
%    experiment where two laser-Doppler velocimeters take measurements
%    at random points along the centreline of the beam.  These measurements
%    consist of a position (x), an incident angle (p) and the magnitude
%    of the velocity along the line of sight (v).
%    The problem is then to fit
% 
%                          2      3                    2     3
%        v = (c + c x + c x  + c x ) cos[ d + d x + d x + d x  - p ]
%              0   1     2      3          0   1     2     3
%            <---- magnitude ----->       <------ phase ----->
% 
%    in the least-squares sense.
% 
%    Source:
%    a modification of an exercize for L. Watson course on LANCELOT in
%    the Spring 1993. Compared to the original proposal, the unnecessary
%    elements were removed as well as an unnecessary constraint on the phase.
% 
%    SIF input: Ph. L. Toint, May 1993, based on a proposal by
%               D. E. Montgomery, Virginia Tech., April 1993.
%    Nonlinear-equations version of VIBRBEAM.SIF, Nick Gould, Jan 2020.
% 
%    classification = 'C-CNOR2-MN-8-30'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'VIBRBEAMNE';

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
        v_('0') = 0;
        v_('1') = 1;
        v_('3') = 3;
        v_('m') = 30;
        v_('x1') = 39.1722;
        v_('x2') = 53.9707;
        v_('x3') = 47.9829;
        v_('x4') = 12.5925;
        v_('x5') = 16.5414;
        v_('x6') = 18.9548;
        v_('x7') = 27.7168;
        v_('x8') = 31.9201;
        v_('x9') = 45.6830;
        v_('x10') = 22.2524;
        v_('x11') = 33.9805;
        v_('x12') = 6.8425;
        v_('x13') = 35.1677;
        v_('x14') = 33.5682;
        v_('x15') = 43.3659;
        v_('x16') = 13.3835;
        v_('x17') = 25.7273;
        v_('x18') = 21.0230;
        v_('x19') = 10.9755;
        v_('x20') = 1.5323;
        v_('x21') = 45.4416;
        v_('x22') = 14.5431;
        v_('x23') = 22.4313;
        v_('x24') = 29.0144;
        v_('x25') = 25.2675;
        v_('x26') = 15.5095;
        v_('x27') = 9.6297;
        v_('x28') = 8.3009;
        v_('x29') = 30.8694;
        v_('x30') = 43.3299;
        v_('v1') = -1.2026;
        v_('v2') = 1.7053;
        v_('v3') = 0.5410;
        v_('v4') = 1.1477;
        v_('v5') = 1.2447;
        v_('v6') = 0.9428;
        v_('v7') = -0.1360;
        v_('v8') = -0.7542;
        v_('v9') = -0.3396;
        v_('v10') = 0.7057;
        v_('v11') = -0.8509;
        v_('v12') = -0.1201;
        v_('v13') = -1.2193;
        v_('v14') = -1.0448;
        v_('v15') = -0.7723;
        v_('v16') = 0.4342;
        v_('v17') = 0.1154;
        v_('v18') = 0.2868;
        v_('v19') = 0.3558;
        v_('v20') = -0.5090;
        v_('v21') = -0.0842;
        v_('v22') = 0.6021;
        v_('v23') = 0.1197;
        v_('v24') = -0.1827;
        v_('v25') = 0.1806;
        v_('v26') = 0.5395;
        v_('v27') = 0.2072;
        v_('v28') = 0.1466;
        v_('v29') = -0.2672;
        v_('v30') = -0.3038;
        v_('p1') = 2.5736;
        v_('p2') = 2.7078;
        v_('p3') = 2.6613;
        v_('p4') = 2.0374;
        v_('p5') = 2.1553;
        v_('p6') = 2.2195;
        v_('p7') = 2.4077;
        v_('p8') = 2.4772;
        v_('p9') = 2.6409;
        v_('p10') = 2.2981;
        v_('p11') = 2.5073;
        v_('p12') = 1.8380;
        v_('p13') = 2.5236;
        v_('p14') = 2.5015;
        v_('p15') = 2.6186;
        v_('p16') = 0.4947;
        v_('p17') = 0.6062;
        v_('p18') = 0.5588;
        v_('p19') = 0.4772;
        v_('p20') = 0.4184;
        v_('p21') = 0.9051;
        v_('p22') = 0.5035;
        v_('p23') = 0.5723;
        v_('p24') = 0.6437;
        v_('p25') = 0.6013;
        v_('p26') = 0.5111;
        v_('p27') = 0.4679;
        v_('p28') = 0.4590;
        v_('p29') = 0.6666;
        v_('p30') = 0.8630;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','c0',ix_);
        pb.xnames{iv} = 'c0';
        [iv,ix_] = s2mpjlib('ii','c1',ix_);
        pb.xnames{iv} = 'c1';
        [iv,ix_] = s2mpjlib('ii','c2',ix_);
        pb.xnames{iv} = 'c2';
        [iv,ix_] = s2mpjlib('ii','c3',ix_);
        pb.xnames{iv} = 'c3';
        [iv,ix_] = s2mpjlib('ii','d0',ix_);
        pb.xnames{iv} = 'd0';
        [iv,ix_] = s2mpjlib('ii','d1',ix_);
        pb.xnames{iv} = 'd1';
        [iv,ix_] = s2mpjlib('ii','d2',ix_);
        pb.xnames{iv} = 'd2';
        [iv,ix_] = s2mpjlib('ii','d3',ix_);
        pb.xnames{iv} = 'd3';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for i=v_('1'):v_('m')
            [ig,ig_] = s2mpjlib('ii',['f',int2str(i)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['f',int2str(i)];
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
        for i=v_('1'):v_('m')
            pbm.gconst(ig_(['f',int2str(i)])) = v_(['v',int2str(i)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_('c0'),1) = -3.5;
        pb.x0(ix_('c1'),1) = 1.0;
        pb.x0(ix_('d0'),1) = 1.7;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'efun',iet_);
        elftv{it}{1} = 'a0';
        elftv{it}{2} = 'a1';
        elftv{it}{3} = 'a2';
        elftv{it}{4} = 'a3';
        elftv{it}{5} = 'b';
        elftp{it}{1} = 'y';
        elftp{it}{2} = 'q';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for i=v_('1'):v_('m')
            for j=v_('0'):v_('3')
                ename = ['fu',int2str(i),',',int2str(j)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'efun';
                ielftype(ie) = iet_('efun');
                vname = 'd0';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('a0',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'd1';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('a1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'd2';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('a2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = 'd3';
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('a3',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['c',int2str(j)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('b',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('y',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['x',int2str(i)]);
                [~,posep] = ismember('q',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['p',int2str(i)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for i=v_('1'):v_('m')
            v_('y') = 1.0;
            for j=v_('0'):v_('3')
                ig = ig_(['f',int2str(i)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['fu',int2str(i),',',int2str(j)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('y');
                v_('y') = v_('y')*v_(['x',int2str(i)]);
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
% LO SOLUTION             0.15644607137
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-8-30';
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

    case 'efun'

        EV_  = varargin{1};
        iel_ = varargin{2};
        y2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        y3 = pbm.elpar{iel_}(1)*y2;
        y4 = y2*y2;
        y5 = y2*y3;
        y6 = y3*y3;
        phi =...
              EV_(1)+pbm.elpar{iel_}(1)*(EV_(2)+pbm.elpar{iel_}(1)*(EV_(3)+pbm.elpar{iel_}(1)*EV_(4)))-pbm.elpar{iel_}(2);
        cosphi = cos(phi);
        sinphi = sin(phi);
        bcos = EV_(5)*cosphi;
        bsin = EV_(5)*sinphi;
        varargout{1} = bcos;
        if(nargout>1)
            g_(1,1) = -bsin;
            g_(2,1) = -bsin*pbm.elpar{iel_}(1);
            g_(3,1) = -bsin*y2;
            g_(4,1) = -bsin*y3;
            g_(5,1) = cosphi;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = -bcos;
                H_(1,2) = -bcos*pbm.elpar{iel_}(1);
                H_(2,1) = H_(1,2);
                H_(1,3) = -bcos*y2;
                H_(3,1) = H_(1,3);
                H_(1,4) = -bcos*y3;
                H_(4,1) = H_(1,4);
                H_(1,5) = -sinphi;
                H_(5,1) = H_(1,5);
                H_(2,2) = -bcos*y2;
                H_(2,3) = -bcos*y3;
                H_(3,2) = H_(2,3);
                H_(2,4) = -bcos*y4;
                H_(4,2) = H_(2,4);
                H_(2,5) = -sinphi*pbm.elpar{iel_}(1);
                H_(5,2) = H_(2,5);
                H_(3,3) = -bcos*y4;
                H_(3,4) = -bcos*y5;
                H_(4,3) = H_(3,4);
                H_(3,5) = -sinphi*y2;
                H_(5,3) = H_(3,5);
                H_(4,4) = -bcos*y6;
                H_(4,5) = -sinphi*y3;
                H_(5,4) = H_(4,5);
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

