function varargout = PDE2(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PDE2
%    *********
% 
%    The pde_2, _20 & _200.mod AMPL models from Hans Mittelmann 
%    (mittelmann@asu.edu)
%    See: http://plato.asu.edu/ftp/barrier/
% 
%    SIF input: Nick Gould, April 25th 2012
%               correction by S. Gratton & Ph. Toint, May 2024
% 
%    classification = 'LLR2-AN-V-V'
% 
%    the x-y discretization 
% 
%       Alternative values for the SIF file parameters:
% IE N                   3              $-PARAMETER
% IE N                   299            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PDE2';

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
            v_('N') = 6;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   2999           $-PARAMETER     pde_2.mod value
% IE N                   2099           $-PARAMETER     pde_20.mod value
% IE N                   1299           $-PARAMETER     pde_200.mod value
        v_('0') = 0;
        v_('1') = 1;
        v_('ONE') = 1.0;
        v_('N1') = 1+v_('N');
        v_('RN1') = v_('N1');
        v_('A') = 0.01;
        v_('G') = 20.0;
        v_('H') = v_('ONE')/v_('RN1');
        v_('-H') = -1.0*v_('H');
        v_('H2') = v_('H')*v_('H');
        v_('GH2') = v_('G')*v_('H2');
        v_('AH') = v_('A')*v_('H');
        v_('SQRTAH') = sqrt(v_('AH'));
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N1')
            for J=v_('0'):v_('N1')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('0'):v_('N1')
                [iv,ix_] = s2mpjlib('ii',['T',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['T',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            for J=v_('0'):v_('N1')
                [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
                gtype{ig} = '<>';
                iv = ix_(['T',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            v_('I+') = 1+I;
            v_('I-') = -1+I;
            for J=v_('1'):v_('N')
                v_('J+') = 1+J;
                v_('J-') = -1+J;
                [ig,ig_] = s2mpjlib('ii',['P',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['P',int2str(I),',',int2str(J)];
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 4.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 4.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('J+')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(round(v_('J-')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('I+'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(round(v_('I-'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['A',int2str(I),',',int2str(J)];
                iv = ix_(['T',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('H');
                end
                [ig,ig_] = s2mpjlib('ii',['B',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '<=';
                cnames{ig} = ['B',int2str(I),',',int2str(J)];
                iv = ix_(['T',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                iv = ix_(['X',int2str(I),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('H')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('H');
                end
            end
        end
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(round(v_('0')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I),',',int2str(round(v_('0')))];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(round(v_('0')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I),',',int2str(round(v_('0')))];
            iv = ix_(['X',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I),',',int2str(round(v_('0')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['D',int2str(I),',',int2str(round(v_('0')))];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I),',',int2str(round(v_('0')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['D',int2str(I),',',int2str(round(v_('0')))];
            iv = ix_(['X',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(round(v_('N1')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I),',',int2str(round(v_('N1')))];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I),',',int2str(round(v_('N1')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I),',',int2str(round(v_('N1')))];
            iv = ix_(['X',int2str(I),',',int2str(round(v_('N1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I),',',int2str(round(v_('N1')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['D',int2str(I),',',int2str(round(v_('N1')))];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['D',int2str(I),',',int2str(round(v_('N1')))],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['D',int2str(I),',',int2str(round(v_('N1')))];
            iv = ix_(['X',int2str(I),',',int2str(round(v_('N1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('0'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['E',int2str(round(v_('0'))),',',int2str(I)];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('0'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['E',int2str(round(v_('0'))),',',int2str(I)];
            iv = ix_(['X',int2str(round(v_('0'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['F',int2str(round(v_('0'))),',',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['F',int2str(round(v_('0'))),',',int2str(I)];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['F',int2str(round(v_('0'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['F',int2str(round(v_('0'))),',',int2str(I)];
            iv = ix_(['X',int2str(round(v_('0'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('N1'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['E',int2str(round(v_('N1'))),',',int2str(I)];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('N1'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['E',int2str(round(v_('N1'))),',',int2str(I)];
            iv = ix_(['X',int2str(round(v_('N1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
            [ig,ig_] = s2mpjlib('ii',['F',int2str(round(v_('N1'))),',',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['F',int2str(round(v_('N1'))),',',int2str(I)];
            iv = ix_(['T',int2str(I),',',int2str(round(v_('0')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            [ig,ig_] = s2mpjlib('ii',['F',int2str(round(v_('N1'))),',',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['F',int2str(round(v_('N1'))),',',int2str(I)];
            iv = ix_(['X',int2str(round(v_('N1'))),',',int2str(I)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('SQRTAH')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('SQRTAH');
            end
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
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('N')
                pbm.gconst(ig_(['P',int2str(I),',',int2str(J)])) = v_('GH2');
            end
        end
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_('IH') = v_('RI')*v_('H');
            v_('IH-1') = -1.0+v_('IH');
            v_('P') = v_('RI')*v_('IH-1');
            v_('P') = 5.0*v_('P');
            for J=v_('1'):v_('N')
                v_('RJ') = J;
                v_('JH') = v_('RJ')*v_('H');
                v_('JH-1') = -1.0+v_('JH');
                v_('YD') = v_('RJ')*v_('JH-1');
                v_('YD') = v_('YD')*v_('P');
                v_('YD') = 3.0+v_('YD');
                v_('YD') = v_('YD')*v_('H');
                v_('-YD') = -1.0*v_('YD');
                pbm.gconst(ig_(['A',int2str(I),',',int2str(J)])) = v_('YD');
                pbm.gconst(ig_(['B',int2str(I),',',int2str(J)])) = v_('-YD');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 3.5*ones(pb.n,1);
        for I=v_('1'):v_('N')
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('0')))])) = 10.0;
            pb.xupper(ix_(['X',int2str(I),',',int2str(round(v_('N1')))])) = 10.0;
            pb.xupper(ix_(['X',int2str(round(v_('0'))),',',int2str(I)])) = 10.0;
            pb.xupper(ix_(['X',int2str(round(v_('N1'))),',',int2str(I)])) = 10.0;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'LLR2-AN-V-V';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
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

