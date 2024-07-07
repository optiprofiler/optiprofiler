function varargout = PALMER7C(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER7C
%    *********
% 
%    A linear least squares problem
%    arising from chemical kinetics.
% 
%     model: H-N=C=Se TZVP + MP2
%    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
%                 A10 X**10 + A12 X**12 + A14 X**14
% 
%    Source:
%    M. Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1992.
% 
%    classification = 'QUR2-RN-8-0'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER7C';

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
        v_('M') = 24;
        v_('1') = 1;
        v_('12') = 12;
        v_('X12') = 0.000000;
        v_('X13') = 0.139626;
        v_('X14') = 0.261799;
        v_('X15') = 0.436332;
        v_('X16') = 0.565245;
        v_('X17') = 0.512942;
        v_('X18') = 0.610865;
        v_('X19') = 0.785398;
        v_('X20') = 0.959931;
        v_('X21') = 1.134464;
        v_('X22') = 1.308997;
        v_('X23') = 1.483530;
        v_('X24') = 1.658063;
        v_('Y12') = 4.419446;
        v_('Y13') = 3.564931;
        v_('Y14') = 2.139067;
        v_('Y15') = 0.404686;
        v_('Y16') = 0.000000;
        v_('Y17') = 0.035152;
        v_('Y18') = 0.146813;
        v_('Y19') = 2.718058;
        v_('Y20') = 9.474417;
        v_('Y21') = 26.132221;
        v_('Y22') = 41.451561;
        v_('Y23') = 72.283164;
        v_('Y24') = 117.630959;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','A0',ix_);
        pb.xnames{iv} = 'A0';
        [iv,ix_] = s2mpjlib('ii','A2',ix_);
        pb.xnames{iv} = 'A2';
        [iv,ix_] = s2mpjlib('ii','A4',ix_);
        pb.xnames{iv} = 'A4';
        [iv,ix_] = s2mpjlib('ii','A6',ix_);
        pb.xnames{iv} = 'A6';
        [iv,ix_] = s2mpjlib('ii','A8',ix_);
        pb.xnames{iv} = 'A8';
        [iv,ix_] = s2mpjlib('ii','A10',ix_);
        pb.xnames{iv} = 'A10';
        [iv,ix_] = s2mpjlib('ii','A12',ix_);
        pb.xnames{iv} = 'A12';
        [iv,ix_] = s2mpjlib('ii','A14',ix_);
        pb.xnames{iv} = 'A14';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('12'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            v_('XQUART') = v_('XSQR')*v_('XSQR');
            v_('X**6') = v_('XSQR')*v_('XQUART');
            v_('X**8') = v_('XSQR')*v_('X**6');
            v_('X**10') = v_('XSQR')*v_('X**8');
            v_('X**12') = v_('XSQR')*v_('X**10');
            v_('X**14') = v_('XSQR')*v_('X**12');
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('A0');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            iv = ix_('A2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XSQR')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XSQR');
            end
            iv = ix_('A4');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XQUART')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XQUART');
            end
            iv = ix_('A6');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**6')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**6');
            end
            iv = ix_('A8');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**8')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**8');
            end
            iv = ix_('A10');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**10')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**10');
            end
            iv = ix_('A12');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**12')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**12');
            end
            iv = ix_('A14');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('X**14')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('X**14');
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('12'):v_('M')
            pbm.gconst(ig_(['O',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('A0')) = -Inf;
        pb.xupper(ix_('A0'),1) = +Inf;
        pb.xlower(ix_('A2')) = -Inf;
        pb.xupper(ix_('A2'),1) = +Inf;
        pb.xlower(ix_('A4')) = -Inf;
        pb.xupper(ix_('A4'),1) = +Inf;
        pb.xlower(ix_('A6')) = -Inf;
        pb.xupper(ix_('A6'),1) = +Inf;
        pb.xlower(ix_('A8')) = -Inf;
        pb.xupper(ix_('A8'),1) = +Inf;
        pb.xlower(ix_('A10')) = -Inf;
        pb.xupper(ix_('A10'),1) = +Inf;
        pb.xlower(ix_('A12')) = -Inf;
        pb.xupper(ix_('A12'),1) = +Inf;
        pb.xlower(ix_('A14')) = -Inf;
        pb.xupper(ix_('A14'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('12'):v_('M')
            ig = ig_(['O',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN              5.0310687D-02
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'QUR2-RN-8-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% ********************
%  SET UP THE GROUPS *
%  ROUTINE           *
% ********************

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

