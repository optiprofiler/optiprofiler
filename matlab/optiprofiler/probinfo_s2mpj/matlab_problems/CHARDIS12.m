function varargout = CHARDIS12(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHARDIS12
%    *********
% 
%    Distribution of charges on a round plate (2D)
% 
%    SIF input: R. Felkel, Jun 1999.
%               correction by S. Gratton & Ph. Toint, May 2024
%    modifield version of CHARDIS1 (formulation corrected)
% 
%    classification = 'OQR2-AY-V-V'
% 
%    Number of positive (or negative) charges -> Number of variables 2*NP1
% 
%       Alternative values for the SIF file parameters:
% IE NP1                 5              $-PARAMETER
% IE NP1                 8              $-PARAMETER
% IE NP1                 20             $-PARAMETER     original value
% IE NP1                 50             $-PARAMETER
% IE NP1                 100            $-PARAMETER
% IE NP1                 200            $-PARAMETER
% IE NP1                 500            $-PARAMETER
% IE NP1                 1000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHARDIS12';

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
            v_('NP1') = 20;  %  SIF file default value
        else
            v_('NP1') = varargin{1};
        end
% IE NP1                 2000           $-PARAMETER
% IE NP1                 5000           $-PARAMETER
        v_('R') = 1.0;
        v_('R2') = v_('R')*v_('R');
        v_('N') = -1+v_('NP1');
        v_('NReal') = v_('N');
        v_('NP1Real') = v_('NP1');
        v_('halfPI') = asin(1.0);
        v_('PI') = 2.0*v_('halfPI');
        v_('2PI') = 4.0*v_('halfPI');
        v_('4PI') = 8.0*v_('halfPI');
        v_('4PIqN') = v_('4PI')/v_('NReal');
        v_('2PIqN') = v_('2PI')/v_('NReal');
        v_('PIqN') = v_('PI')/v_('NReal');
        v_('RqN') = v_('R')/v_('NReal');
        v_('1') = 1;
        v_('2') = 2;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('NP1')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Y',int2str(I)],ix_);
            pb.xnames{iv} = ['Y',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('NP1')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP1')
                [ig,ig_] = s2mpjlib('ii',['O',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        for I=v_('2'):v_('NP1')
            [ig,ig_] = s2mpjlib('ii',['RES',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['RES',int2str(I)];
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
        for I=v_('2'):v_('NP1')
            pbm.gconst(ig_(['RES',int2str(I)])) = v_('R2');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('2'):v_('NP1')
            pb.xlower(ix_(['X',int2str(I)])) = -Inf;
            pb.xupper(ix_(['X',int2str(I)]),1) = +Inf;
            pb.xlower(ix_(['Y',int2str(I)])) = -Inf;
            pb.xupper(ix_(['Y',int2str(I)]),1) = +Inf;
        end
        pb.xlower(ix_('X1'),1) = v_('R');
        pb.xupper(ix_('X1'),1) = v_('R');
        pb.xlower(ix_('Y1'),1) = 0.0;
        pb.xupper(ix_('Y1'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('2'):v_('NP1')
            v_('I-') = -1+I;
            v_('RealI-') = v_('I-');
            v_('RealNP1-I') = v_('NP1Real')-v_('RealI-');
            v_('PHII-') = v_('2PIqN')*v_('RealI-');
            v_('RI-') = v_('RqN')*v_('RealNP1-I');
            v_('XST') = cos(v_('PHII-'));
            v_('YST') = sin(v_('PHII-'));
            v_('XS') = v_('XST')*v_('RI-');
            v_('YS') = v_('YST')*v_('RI-');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('XS');
            pb.x0(ix_(['Y',int2str(I)]),1) = v_('YS');
        end
        pb.x0(ix_('X1'),1) = v_('R');
        pb.x0(ix_('Y1'),1) = 0.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'V1';
        [it,iet_] = s2mpjlib( 'ii', 'eDIFSQR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NP1')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP1')
                ename = ['X',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eDIFSQR';
                ielftype(ie) = iet_('eDIFSQR');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['Y',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eDIFSQR';
                ielftype(ie) = iet_('eDIFSQR');
                vname = ['Y',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['Y',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                posev = find(strcmp('V2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        for I=v_('2'):v_('NP1')
            ename = ['RX',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['RY',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['Y',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gREZIP',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NP1')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP1')
                ig = ig_(['O',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gREZIP';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['X',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['Y',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        for I=v_('2'):v_('NP1')
            ig = ig_(['RES',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['RX',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['RY',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OQR2-AY-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

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

    case 'eDIFSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-EV_(2))*(EV_(1)-EV_(2));
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-EV_(2));
            g_(2,1) = -2.0*(EV_(1)-EV_(2));
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0;
                H_(1,2) = -2.0;
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gREZIP'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = 1.0/GVAR_;
        if(nargout>1)
            g_ = -1.0/(GVAR_*GVAR_);
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0/(GVAR_*GVAR_*GVAR_);
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

