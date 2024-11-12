function varargout = CHARDIS02(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : CHARDIS02
%    *********
% 
%    Distribution of (equal)charges on [-R,R]x[-R,R] (2D)
% 
%    SIF input: R. Felkel, Jun 1999.
%               correction by S. Gratton & Ph. Toint, May 2024
%    modifield version of CHARDIS0 (formulation corrected)
% 
%    classification = 'C-COBR2-AY-V-V'
% 
%    Number of positive (or negative) charges -> Number of variables 2*NP1
% 
%       Alternative values for the SIF file parameters:
% IE NP1                 5              $-PARAMETER
% IE NP1                 9              $-PARAMETER
% IE NP1                 20             $-PARAMETER
% IE NP1                 30             $-PARAMETER
% IE NP1                 50             $-PARAMETER     original value
% IE NP1                 100            $-PARAMETER
% IE NP1                 200            $-PARAMETER
% IE NP1                 500            $-PARAMETER
% IE NP1                 1000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'CHARDIS02';

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
        if(nargs<1)
            v_('NP1') = 20;  %  SIF file default value
        else
            v_('NP1') = varargin{1};
        end
% IE NP1                 2000           $-PARAMETER
% IE NP1                 5000           $-PARAMETER
        v_('R') = 10.0;
        v_('R-') = -10.0;
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
                pbm.gscale(ig,1) = 0.01;
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('1'):v_('NP1')
            pb.xlower(ix_(['X',int2str(I)]),1) = v_('R-');
            pb.xupper(ix_(['X',int2str(I)])) = v_('R');
            pb.xlower(ix_(['Y',int2str(I)]),1) = v_('R-');
            pb.xupper(ix_(['Y',int2str(I)])) = v_('R');
        end
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        for I=v_('1'):v_('NP1')
            v_('RealI-') = I;
            v_('RealNP1-I') = v_('NP1Real')-v_('RealI-');
            v_('PHII-') = v_('2PIqN')*v_('RealI-');
            v_('RI-') = v_('RqN')*v_('RealNP1-I');
            v_('XSTT') = cos(v_('PHII-'));
            v_('YSTT') = sin(v_('PHII-'));
            v_('XST') = v_('XSTT')*v_('RI-');
            v_('YST') = v_('YSTT')*v_('RI-');
            v_('XS') = 0.5*v_('XST');
            v_('YS') = 0.5*v_('YST');
            pb.x0(ix_(['X',int2str(I)]),1) = v_('XS');
            pb.x0(ix_(['Y',int2str(I)]),1) = v_('YS');
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eDIFSQR',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
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
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gREZIP',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('NP1')
            v_('I+') = 1+I;
            for J=v_('I+'):v_('NP1')
                ig = ig_(['O',int2str(I),',',int2str(J)]);
                pbm.grftype{ig} = 'gREZIP';
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['X',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['Y',int2str(I),',',int2str(J)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-COBR2-AY-V-V';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

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

