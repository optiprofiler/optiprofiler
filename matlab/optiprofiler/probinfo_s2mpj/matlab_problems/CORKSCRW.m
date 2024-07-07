function varargout = CORKSCRW(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CORKSCRW
%    *********
% 
%    A nonlinear optimal control problem with both state- and
%    control constraints.
%    The problem is to control (using an applied force of limited
%    magnitude) a mass moving in the 3D space, such that its
%    trajectory lies within a prescribed distance TOL of the
%    corkscreww-like curve defined by
%               y = sin(x), z = cos(x),
%    and such that it stops at a given abscissa XT in minimum time.
%    The mass is initially stationary at (0,0,1).
% 
%    Source:
%    Ph. Toint, private communication.
% 
%    SIF input: Ph. Toint, April 1991.
% 
%    classification = 'SOR2-AN-V-V'
% 
%    Number of time intervals
%    The number of variables is 9T+6, of which 9 are fixed.
% 
%       Alternative values for the SIF file parameters:
% IE T                   10             $-PARAMETER n = 96     original value
% IE T                   50             $-PARAMETER n = 456
% IE T                   100            $-PARAMETER n = 906
% IE T                   500            $-PARAMETER n = 4506
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CORKSCRW';

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
            v_('T') = 10;  %  SIF file default value
        else
            v_('T') = varargin{1};
        end
% IE T                   1000           $-PARAMETER n = 9006
        if(nargs<2)
            v_('XT') = 10.0;  %  SIF file default value
        else
            v_('XT') = varargin{2};
        end
        if(nargs<3)
            v_('MASS') = 0.37;  %  SIF file default value
        else
            v_('MASS') = varargin{3};
        end
        if(nargs<4)
            v_('TOL') = 0.1;  %  SIF file default value
        else
            v_('TOL') = varargin{4};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('RT') = v_('T');
        v_('T+1') = 1.0+v_('RT');
        v_('H') = v_('XT')/v_('RT');
        v_('1/H') = 1.0/v_('H');
        v_('-1/H') = -1.0*v_('1/H');
        v_('M/H') = v_('MASS')/v_('H');
        v_('-M/H') = -1.0*v_('M/H');
        v_('TOLSQ') = v_('TOL')*v_('TOL');
        v_('XTT+1') = v_('XT')*v_('T+1');
        v_('W') = 0.5*v_('XTT+1');
        for I=v_('1'):v_('T')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            v_(['W/T',int2str(I)]) = v_('W')/v_('TI');
        end
        v_('FMAX') = v_('XT')/v_('RT');
        v_('-FMAX') = -1.0*v_('FMAX');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('T')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Z',int2str(I)],ix_);
            pb.xnames{iv} = ['Z',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['VX',int2str(I)],ix_);
            pb.xnames{iv} = ['VX',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['VY',int2str(I)],ix_);
            pb.xnames{iv} = ['VY',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['VZ',int2str(I)],ix_);
            pb.xnames{iv} = ['VZ',int2str(I)];
        end
        for I=v_('1'):v_('T')
            [iv,ix_] = s2mpjlib('ii',['UX',int2str(I)],ix_);
            pb.xnames{iv} = ['UX',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['UY',int2str(I)],ix_);
            pb.xnames{iv} = ['UY',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['UZ',int2str(I)],ix_);
            pb.xnames{iv} = ['UZ',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('T')
            [ig,ig_] = s2mpjlib('ii',['OX',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = v_(['W/T',int2str(I)]);
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
        end
        for I=v_('1'):v_('T')
            v_('I-1') = -1+I;
            [ig,ig_] = s2mpjlib('ii',['ACX',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ACX',int2str(I)];
            iv = ix_(['VX',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('M/H');
            end
            iv = ix_(['VX',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-M/H');
            end
            iv = ix_(['UX',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ACY',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ACY',int2str(I)];
            iv = ix_(['VY',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('M/H');
            end
            iv = ix_(['VY',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-M/H');
            end
            iv = ix_(['UY',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['ACZ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['ACZ',int2str(I)];
            iv = ix_(['VZ',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('M/H');
            end
            iv = ix_(['VZ',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-M/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-M/H');
            end
            iv = ix_(['UZ',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['PSX',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PSX',int2str(I)];
            iv = ix_(['X',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['X',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['VX',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['PSY',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PSY',int2str(I)];
            iv = ix_(['Y',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['Y',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['VY',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['PSZ',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PSZ',int2str(I)];
            iv = ix_(['Z',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('1/H');
            end
            iv = ix_(['Z',int2str(round(v_('I-1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-1/H')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-1/H');
            end
            iv = ix_(['VZ',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['SC',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['SC',int2str(I)];
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
        for I=v_('1'):v_('T')
            pbm.gconst(ig_(['OX',int2str(I)])) = v_('XT');
            pbm.gconst(ig_(['SC',int2str(I)])) = v_('TOLSQ');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        pb.xlower(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['Y',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['Y',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['Z',int2str(round(v_('0')))]),1) = 1.0;
        pb.xupper(ix_(['Z',int2str(round(v_('0')))]),1) = 1.0;
        pb.xlower(ix_(['VX',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['VX',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['VY',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['VY',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['VZ',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['VZ',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['VX',int2str(round(v_('T')))]),1) = 0.0;
        pb.xupper(ix_(['VX',int2str(round(v_('T')))]),1) = 0.0;
        pb.xlower(ix_(['VY',int2str(round(v_('T')))]),1) = 0.0;
        pb.xupper(ix_(['VY',int2str(round(v_('T')))]),1) = 0.0;
        pb.xlower(ix_(['VZ',int2str(round(v_('T')))]),1) = 0.0;
        pb.xupper(ix_(['VZ',int2str(round(v_('T')))]),1) = 0.0;
        for I=v_('1'):v_('T')
            pb.xlower(ix_(['UX',int2str(I)]),1) = v_('-FMAX');
            pb.xupper(ix_(['UX',int2str(I)])) = v_('FMAX');
            pb.xlower(ix_(['UY',int2str(I)]),1) = v_('-FMAX');
            pb.xupper(ix_(['UY',int2str(I)])) = v_('FMAX');
            pb.xlower(ix_(['UZ',int2str(I)]),1) = v_('-FMAX');
            pb.xupper(ix_(['UZ',int2str(I)])) = v_('FMAX');
            pb.xlower(ix_(['X',int2str(I)]),1) = 0.0;
            pb.xupper(ix_(['X',int2str(I)])) = v_('XT');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        pb.x0(ix_(['X',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['Y',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['Z',int2str(round(v_('0')))]),1) = 1.0;
        pb.x0(ix_(['VX',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['VY',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['VZ',int2str(round(v_('0')))]),1) = 0.0;
        pb.x0(ix_(['VX',int2str(round(v_('T')))]),1) = 0.0;
        pb.x0(ix_(['VY',int2str(round(v_('T')))]),1) = 0.0;
        pb.x0(ix_(['VZ',int2str(round(v_('T')))]),1) = 0.0;
        for I=v_('1'):v_('T')
            v_('RI') = I;
            v_('TI') = v_('RI')*v_('H');
            if(isKey(ix_,['X',int2str(I)]))
                pb.x0(ix_(['X',int2str(I)]),1) = v_('TI');
            else
                pb.y0(find(pbm.congrps==ig_(['X',int2str(I)])),1) = v_('TI');
            end
            if(isKey(ix_,['VX',int2str(I)]))
                pb.x0(ix_(['VX',int2str(I)]),1) = 1.0;
            else
                pb.y0(find(pbm.congrps==ig_(['VX',int2str(I)])),1) = 1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eERRSIN',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eERRCOS',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Z';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('T')
            ename = ['ES',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eERRSIN';
            ielftype(ie) = iet_('eERRSIN');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eERRCOS';
            ielftype(ie) = iet_('eERRCOS');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['Z',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Z',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('T')
            ig = ig_(['OX',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['SC',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ES',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EC',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           1.1601050195
% LO SOLTN(50)           26.484181830
% LO SOLTN(100)          44.368110588
% LO SOLTN(500)
% LO SOLTN(1000)
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'SOR2-AN-V-V';
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

    case 'eERRSIN'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SINX = sin(EV_(1));
        COSX = cos(EV_(1));
        ERR = EV_(2)-SINX;
        varargout{1} = ERR*ERR;
        if(nargout>1)
            g_(1,1) = -2.0*ERR*COSX;
            g_(2,1) = 2.0*ERR;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*(COSX^2+ERR*SINX);
                H_(1,2) = -2.0*COSX;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eERRCOS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        SINX = sin(EV_(1));
        COSX = cos(EV_(1));
        ERR = EV_(2)-COSX;
        varargout{1} = ERR*ERR;
        if(nargout>1)
            g_(1,1) = 2.0*ERR*SINX;
            g_(2,1) = 2.0*ERR;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*(SINX^2+ERR*COSX);
                H_(1,2) = 2.0*SINX;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0;
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

