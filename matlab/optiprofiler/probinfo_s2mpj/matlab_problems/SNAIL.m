function varargout = SNAIL(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SNAIL
%    *********
% 
%    A 2D problem featuring a spiraling valley.
%    Dedicated to the city of Namur, whose emblem is a snail.
% 
%    Source:
%    J. Engels, private communication.
% 
%    SIF input: Ph. Toint, May 1990.
% 
%    classification = 'OUR2-AN-2-0'
% 
%    Problem parameters (CUP > CLOW > 0)
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SNAIL';

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
        v_('CLOW') = 1.0;
        v_('CUP') = 2.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 10.0;
        pb.x0(ix_('X2'),1) = 10.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSPIRAL',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'CL';
        elftp{it}{2} = 'CU';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        ename = 'E';
        [ie,ie_] = s2mpjlib('ii',ename,ie_);
        pbm.elftype{ie} = 'eSPIRAL';
        ielftype(ie) = iet_('eSPIRAL');
        vname = 'X1';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('X',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        vname = 'X2';
        [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
        posev = find(strcmp('Y',elftv{ielftype(ie)}));
        pbm.elvar{ie}(posev) = iv;
        [~,posep] = ismember('CL',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('CLOW');
        [~,posep] = ismember('CU',elftp{ielftype(ie)});
        pbm.elpar{ie}(posep) = v_('CUP');
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('E');
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OUR2-AN-2-0';
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

    case 'eSPIRAL'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = 0.5*(pbm.elpar{iel_}(2)+pbm.elpar{iel_}(1));
        B = 0.5*(pbm.elpar{iel_}(2)-pbm.elpar{iel_}(1));
        X2 = EV_(1)*EV_(1);
        Y2 = EV_(2)*EV_(2);
        R2 = X2+Y2;
        D = 1.0+R2;
        D2 = D*D;
        D3 = D2*D;
        U = R2/D;
        DUDX = (EV_(1)+EV_(1))/D2;
        DUDY = (EV_(2)+EV_(2))/D2;
        D2UDX2 = 2.0*(D-4.0*X2)/D3;
        D2UDY2 = 2.0*(D-4.0*Y2)/D3;
        D2UDXY = -8.0*EV_(1)*EV_(2)/D3;
        THETA = atan2(EV_(2),EV_(1));
        DTDX = -EV_(2)/R2;
        DTDY = EV_(1)/R2;
        R4 = R2*R2;
        D2TDX2 = 2.0*EV_(1)*EV_(2)/R4;
        D2TDY2 = -2.0*EV_(2)*EV_(1)/R4;
        D2TDXY = (Y2-X2)/R4;
        R = sqrt(R2);
        R3 = R*R2;
        DRDX = EV_(1)/R;
        DRDY = EV_(2)/R;
        D2RDX2 = Y2/R3;
        D2RDY2 = X2/R3;
        D2RDXY = -EV_(1)*EV_(2)/R3;
        ARG = R-THETA;
        S = B*sin(ARG);
        C = B*cos(ARG);
        DCDX = -S*(DRDX-DTDX);
        DCDY = -S*(DRDY-DTDY);
        D2CDX2 = -C*(DRDX-DTDX)^2-S*(D2RDX2-D2TDX2);
        D2CDY2 = -C*(DRDY-DTDY)^2-S*(D2RDY2-D2TDY2);
        D2CDXY = -C*(DRDX-DTDX)*(DRDY-DTDY)-S*(D2RDXY-D2TDXY);
        V = 1.0+A*R-R*C;
        DVDX = A*DRDX-DRDX*C-R*DCDX;
        DVDY = A*DRDY-DRDY*C-R*DCDY;
        D2VDX2 = A*D2RDX2-D2RDX2*C-2.0*DRDX*DCDX-R*D2CDX2;
        D2VDY2 = A*D2RDY2-D2RDY2*C-2.0*DRDY*DCDY-R*D2CDY2;
        D2VDXY = A*D2RDXY-D2RDXY*C-DRDX*DCDY-DRDY*DCDX-R*D2CDXY;
        varargout{1} = U*V;
        if(nargout>1)
            g_(1,1) = DUDX*V+U*DVDX;
            g_(2,1) = DUDY*V+U*DVDY;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = D2UDX2*V+2.0*DUDX*DVDX+U*D2VDX2;
                H_(1,2) = D2UDXY*V+DUDX*DVDY+DUDY*DVDX+U*D2VDXY;
                H_(2,1) = H_(1,2);
                H_(2,2) = D2UDY2*V+2.0*DUDY*DVDY+U*D2VDY2;
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

