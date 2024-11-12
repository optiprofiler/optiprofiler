function varargout = HS119(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HS119
%    *********
% 
%    Source: problem 119 in
%    W. Hock and K. Schittkowski,
%    "Test examples for nonlinear programming codes",
%    Lectures Notes in Economics and Mathematical Systems 187, Springer
%    Verlag, Heidelberg, 1981.
% 
%    Original Source: problem 7 in
%    A.R. Colville
%    "A comparative study on nonlinear programming"
%    IBM Scientific Center Report 320-2949, New York, 1968.
% 
%    SIF input: A.R. Conn, March 1991.
% 
%    classification = 'C-COLR2-AN-16-8'
% 
%    Set useful parameters
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HS119';

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
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('9') = 9;
        v_('10') = 10;
        v_('11') = 11;
        v_('12') = 12;
        v_('13') = 13;
        v_('14') = 14;
        v_('15') = 15;
        v_('16') = 16;
        for I=v_('1'):v_('16')
            for J=v_('1'):v_('16')
                v_(['A',int2str(I),',',int2str(J)]) = 0.0;
            end
        end
        for I=v_('1'):v_('8')
            for J=v_('1'):v_('16')
                v_(['B',int2str(I),',',int2str(J)]) = 0.0;
            end
        end
        for I=v_('1'):v_('16')
            v_(['A',int2str(I),',',int2str(I)]) = 1.0;
        end
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('4')))]) = 1.0;
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('7')))]) = 1.0;
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('8')))]) = 1.0;
        v_(['A',int2str(round(v_('1'))),',',int2str(round(v_('16')))]) = 1.0;
        v_(['A',int2str(round(v_('2'))),',',int2str(round(v_('3')))]) = 1.0;
        v_(['A',int2str(round(v_('2'))),',',int2str(round(v_('7')))]) = 1.0;
        v_(['A',int2str(round(v_('2'))),',',int2str(round(v_('10')))]) = 1.0;
        v_(['A',int2str(round(v_('3'))),',',int2str(round(v_('7')))]) = 1.0;
        v_(['A',int2str(round(v_('3'))),',',int2str(round(v_('9')))]) = 1.0;
        v_(['A',int2str(round(v_('3'))),',',int2str(round(v_('10')))]) = 1.0;
        v_(['A',int2str(round(v_('3'))),',',int2str(round(v_('14')))]) = 1.0;
        v_(['A',int2str(round(v_('4'))),',',int2str(round(v_('7')))]) = 1.0;
        v_(['A',int2str(round(v_('4'))),',',int2str(round(v_('11')))]) = 1.0;
        v_(['A',int2str(round(v_('4'))),',',int2str(round(v_('15')))]) = 1.0;
        v_(['A',int2str(round(v_('5'))),',',int2str(round(v_('6')))]) = 1.0;
        v_(['A',int2str(round(v_('5'))),',',int2str(round(v_('10')))]) = 1.0;
        v_(['A',int2str(round(v_('5'))),',',int2str(round(v_('12')))]) = 1.0;
        v_(['A',int2str(round(v_('5'))),',',int2str(round(v_('16')))]) = 1.0;
        v_(['A',int2str(round(v_('6'))),',',int2str(round(v_('8')))]) = 1.0;
        v_(['A',int2str(round(v_('6'))),',',int2str(round(v_('15')))]) = 1.0;
        v_(['A',int2str(round(v_('7'))),',',int2str(round(v_('11')))]) = 1.0;
        v_(['A',int2str(round(v_('7'))),',',int2str(round(v_('13')))]) = 1.0;
        v_(['A',int2str(round(v_('8'))),',',int2str(round(v_('10')))]) = 1.0;
        v_(['A',int2str(round(v_('8'))),',',int2str(round(v_('15')))]) = 1.0;
        v_(['A',int2str(round(v_('9'))),',',int2str(round(v_('12')))]) = 1.0;
        v_(['A',int2str(round(v_('9'))),',',int2str(round(v_('16')))]) = 1.0;
        v_(['A',int2str(round(v_('10'))),',',int2str(round(v_('14')))]) = 1.0;
        v_(['A',int2str(round(v_('11'))),',',int2str(round(v_('13')))]) = 1.0;
        v_(['A',int2str(round(v_('12'))),',',int2str(round(v_('14')))]) = 1.0;
        v_(['A',int2str(round(v_('13'))),',',int2str(round(v_('14')))]) = 1.0;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('1')))]) = 0.22;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('1')))]) = -1.46;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('1')))]) = 1.29;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('1')))]) = -1.10;
        v_(['B',int2str(round(v_('7'))),',',int2str(round(v_('1')))]) = 1.12;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('2')))]) = 0.20;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('2')))]) = -0.89;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('2')))]) = -1.06;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('2')))]) = -1.72;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('2')))]) = 0.45;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('3')))]) = 0.19;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('3')))]) = -1.30;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('3')))]) = 0.95;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('3')))]) = -0.33;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('3')))]) = 0.26;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('4')))]) = 0.25;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('4')))]) = 1.82;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('4')))]) = -0.54;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('4')))]) = -1.43;
        v_(['B',int2str(round(v_('7'))),',',int2str(round(v_('4')))]) = 0.31;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('4')))]) = -1.10;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('5')))]) = 0.15;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('5')))]) = -1.15;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('5')))]) = -1.16;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('5')))]) = 1.51;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('5')))]) = 1.62;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('5')))]) = 0.58;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('6')))]) = 0.11;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('6')))]) = -0.96;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('6')))]) = -1.78;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('6')))]) = 0.59;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('6')))]) = 1.24;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('7')))]) = 0.12;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('7')))]) = 0.80;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('7')))]) = -0.41;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('7')))]) = -0.33;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('7')))]) = 0.21;
        v_(['B',int2str(round(v_('7'))),',',int2str(round(v_('7')))]) = 1.12;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('7')))]) = -1.03;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('8')))]) = 0.13;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('8')))]) = -0.49;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('8')))]) = -0.43;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('8')))]) = -0.26;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('8')))]) = 0.10;
        v_(['B',int2str(round(v_('1'))),',',int2str(round(v_('9')))]) = 1.00;
        v_(['B',int2str(round(v_('7'))),',',int2str(round(v_('9')))]) = -0.36;
        v_(['B',int2str(round(v_('2'))),',',int2str(round(v_('10')))]) = 1.00;
        v_(['B',int2str(round(v_('3'))),',',int2str(round(v_('11')))]) = 1.00;
        v_(['B',int2str(round(v_('4'))),',',int2str(round(v_('12')))]) = 1.00;
        v_(['B',int2str(round(v_('5'))),',',int2str(round(v_('13')))]) = 1.00;
        v_(['B',int2str(round(v_('6'))),',',int2str(round(v_('14')))]) = 1.00;
        v_(['B',int2str(round(v_('7'))),',',int2str(round(v_('15')))]) = 1.00;
        v_(['B',int2str(round(v_('8'))),',',int2str(round(v_('16')))]) = 1.00;
        v_(['C',int2str(round(v_('1')))]) = 2.5;
        v_(['C',int2str(round(v_('2')))]) = 1.1;
        v_(['C',int2str(round(v_('3')))]) = -3.1;
        v_(['C',int2str(round(v_('4')))]) = -3.5;
        v_(['C',int2str(round(v_('5')))]) = 1.3;
        v_(['C',int2str(round(v_('6')))]) = 2.1;
        v_(['C',int2str(round(v_('7')))]) = 2.3;
        v_(['C',int2str(round(v_('8')))]) = -1.5;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('16')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('16')
            [ig,ig_] = s2mpjlib('ii',['OG',int2str(I)],ig_);
            gtype{ig} = '<>';
        end
        for I=v_('1'):v_('8')
            for J=v_('1'):v_('16')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['G',int2str(I)];
                iv = ix_(['X',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['B',int2str(I),',',int2str(J)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['B',int2str(I),',',int2str(J)]);
                end
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
        for I=v_('1'):v_('8')
            pbm.gconst(ig_(['G',int2str(I)])) = v_(['C',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xupper = 5.0*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 10.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'U1';
        elftv{it}{2} = 'U2';
        elftp{it}{1} = 'AIJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('16')
            for J=v_('1'):v_('16')
                ename = ['S',int2str(I),',',int2str(J)];
                [ie,ie_,newelt] = s2mpjlib('ii',ename,ie_);
                if(newelt)
                    pbm.elftype{ie} = 'ePROD';
                    ielftype(ie) = iet_('ePROD');
                end
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],5.0,10.0);
                posev = find(strcmp('U1',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],5.0,10.0);
                posev = find(strcmp('U2',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                [~,posep] = ismember('AIJ',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_(['A',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('16')
            for J=v_('1'):v_('16')
                ig = ig_(['OG',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(I),',',int2str(J)]);
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
        pb.pbclass = 'C-COLR2-AN-16-8';
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

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TU1P1 = 2.0*EV_(1)+1;
        TU2P1 = 2.0*EV_(2)+1;
        FIRST = EV_(1)^2+EV_(1)+1.0;
        SECOND = EV_(2)^2+EV_(2)+1.0;
        varargout{1} = pbm.elpar{iel_}(1)*FIRST*SECOND;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*TU1P1*SECOND;
            g_(2,1) = pbm.elpar{iel_}(1)*TU2P1*FIRST;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.elpar{iel_}(1)*2.0*SECOND;
                H_(1,2) = pbm.elpar{iel_}(1)*TU1P1*TU2P1;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*2.0*FIRST;
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

