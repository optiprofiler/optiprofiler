function varargout = DECONVBNE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DECONVBNE
%    *********
% 
%    A problem arising in deconvolution analysis
%    (bounded variables version).
% 
%    Source:
%    J.P. Rasson, Private communication, 1996.
% 
%    SIF input: Ph. Toint, Nov 1996.
%    unititialized variables fixed at zero, Nick Gould, Feb, 2013
%    Bound-constrained nonlinear equations version: Nick Gould, June 2019.
% 
%    classification = 'C-CNOR2-MN-61-0'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DECONVBNE';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('0') = 0;
        v_('1') = 1;
        v_('LGSG') = 11;
        v_('LGTR') = 40;
        v_('-LGSG') = -1*v_('LGSG');
        v_('PIC') = 3.0000000000;
        v_('TR1') = 0.0000000000;
        v_('TR2') = 0.0000000000;
        v_('TR3') = 1.600000E-03;
        v_('TR4') = 5.400000E-03;
        v_('TR5') = 7.020000E-02;
        v_('TR6') = 0.1876000000;
        v_('TR7') = 0.3320000000;
        v_('TR8') = 0.7640000000;
        v_('TR9') = 0.9320000000;
        v_('TR10') = 0.8120000000;
        v_('TR11') = 0.3464000000;
        v_('TR12') = 0.2064000000;
        v_('TR13') = 8.300000E-02;
        v_('TR14') = 3.400000E-02;
        v_('TR15') = 6.179999E-02;
        v_('TR16') = 1.2000000000;
        v_('TR17') = 1.8000000000;
        v_('TR18') = 2.4000000000;
        v_('TR19') = 9.0000000000;
        v_('TR20') = 2.4000000000;
        v_('TR21') = 1.8010000000;
        v_('TR22') = 1.3250000000;
        v_('TR23') = 7.620000E-02;
        v_('TR24') = 0.2104000000;
        v_('TR25') = 0.2680000000;
        v_('TR26') = 0.5520000000;
        v_('TR27') = 0.9960000000;
        v_('TR28') = 0.3600000000;
        v_('TR29') = 0.2400000000;
        v_('TR30') = 0.1510000000;
        v_('TR31') = 2.480000E-02;
        v_('TR32') = 0.2432000000;
        v_('TR33') = 0.3602000000;
        v_('TR34') = 0.4800000000;
        v_('TR35') = 1.8000000000;
        v_('TR36') = 0.4800000000;
        v_('TR37') = 0.3600000000;
        v_('TR38') = 0.2640000000;
        v_('TR39') = 6.000000E-03;
        v_('TR40') = 6.000000E-03;
        v_('SSG1') = 1.000000E-02;
        v_('SSG2') = 2.000000E-02;
        v_('SSG3') = 0.4000000000;
        v_('SSG4') = 0.6000000000;
        v_('SSG5') = 0.8000000000;
        v_('SSG6') = 3.0000000000;
        v_('SSG7') = 0.8000000000;
        v_('SSG8') = 0.6000000000;
        v_('SSG9') = 0.4400000000;
        v_('SSG10') = 1.000000E-02;
        v_('SSG11') = 1.000000E-02;
        v_('CC1') = 0.0;
        v_('CC2') = 0.0;
        v_('CC3') = 0.0;
        v_('CC4') = 0.0;
        v_('CC5') = 0.0;
        v_('CC6') = 0.0;
        v_('CC7') = 0.0;
        v_('CC8') = 0.0;
        v_('CC9') = 0.0;
        v_('CC10') = 0.0;
        v_('CC11') = 0.0;
        v_('CC12') = 0.0;
        v_('CC13') = 0.0;
        v_('CC14') = 0.0;
        v_('CC15') = 0.0;
        v_('CC16') = 0.0;
        v_('CC17') = 0.0;
        v_('CC18') = 0.0;
        v_('CC19') = 0.0;
        v_('CC20') = 0.0;
        v_('CC21') = 0.0;
        v_('CC22') = 0.0;
        v_('CC23') = 0.0;
        v_('CC24') = 0.0;
        v_('CC25') = 0.0;
        v_('CC26') = 0.0;
        v_('CC27') = 0.0;
        v_('CC28') = 0.0;
        v_('CC29') = 0.0;
        v_('CC30') = 0.0;
        v_('CC31') = 0.0;
        v_('CC32') = 0.0;
        v_('CC33') = 0.0;
        v_('CC34') = 0.0;
        v_('CC35') = 0.0;
        v_('CC36') = 0.0;
        v_('CC37') = 0.0;
        v_('CC38') = 0.0;
        v_('CC39') = 0.0;
        v_('CC40') = 0.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for K=v_('-LGSG'):v_('LGTR')
            [iv,ix_] = s2mpjlib('ii',['C',int2str(K)],ix_);
            pb.xnames{iv} = ['C',int2str(K)];
        end
        for I=v_('1'):v_('LGSG')
            [iv,ix_] = s2mpjlib('ii',['SG',int2str(I)],ix_);
            pb.xnames{iv} = ['SG',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for K=v_('1'):v_('LGTR')
            [ig,ig_] = s2mpjlib('ii',['R',int2str(K)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(K)];
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
        for K=v_('1'):v_('LGTR')
            pbm.gconst(ig_(['R',int2str(K)])) = v_(['TR',int2str(K)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('LGSG')
            pb.xlower(ix_(['SG',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['SG',int2str(I)])) = v_('PIC');
        end
        for K=v_('-LGSG'):v_('0')
            pb.xlower(ix_(['C',int2str(K)]),1) = 0.0;
            pb.xupper(ix_(['C',int2str(K)]),1) = 0.0;
        end
        for K=v_('1'):v_('LGTR')
            pb.xlower(ix_(['C',int2str(K)]),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for K=v_('1'):v_('LGTR')
            pb.x0(ix_(['C',int2str(K)]),1) = v_(['CC',int2str(K)]);
        end
        for I=v_('1'):v_('LGSG')
            pb.x0(ix_(['SG',int2str(I)]),1) = v_(['SSG',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'IDX';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for K=v_('1'):v_('LGTR')
            for I=v_('1'):v_('LGSG')
                v_('K-I') = K-I;
                v_('K-I+1') = 1+v_('K-I');
                v_('RIDX') = v_('K-I+1');
                ename = ['PROD',int2str(K),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'ePR';
                ielftype(ie) = iet_('ePR');
                vname = ['SG',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['C',int2str(round(v_('K-I+1')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('IDX',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('RIDX');
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for K=v_('1'):v_('LGTR')
            for I=v_('1'):v_('LGSG')
                ig = ig_(['R',int2str(K)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['PROD',int2str(K),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-61-0';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'ePR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        NEGIDX = pbm.elpar{iel_}(1)<=0.0;
        if(NEGIDX)
            SCAL = 0.0;
        end
        if(~NEGIDX)
            SCAL = 1.0;
        end
        varargout{1} = SCAL*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = SCAL*EV_(2);
            g_(2,1) = SCAL*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = SCAL;
                H_(2,1) = H_(1,2);
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

