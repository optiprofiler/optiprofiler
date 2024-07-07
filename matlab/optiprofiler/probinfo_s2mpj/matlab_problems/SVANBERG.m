function varargout = SVANBERG(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SVANBERG
%    *********
% 
%    A subproblem of the type arising in structural optimization
% 
%    Source:
%    Svanberg K.,
%    "Method of moving asymptots - a new method for structural optimization",
%    Int.J. Num. Meth. Eng, 24, pp. 359--373, 1987
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'OOR2-MN-V-V'
% 
%    Number of variables (must be even and >= 10)
% 
%       Alternative values for the SIF file parameters:
% IE N                   10             $-PARAMETER     original value
% IE N                   20             $-PARAMETER
% IE N                   30             $-PARAMETER
% IE N                   40             $-PARAMETER
% IE N                   50             $-PARAMETER
% IE N                   60             $-PARAMETER
% IE N                   70             $-PARAMETER
% IE N                   80             $-PARAMETER
% IE N                   90             $-PARAMETER
% IE N                   100            $-PARAMETER
% IE N                   500            $-PARAMETER
% IE N                   1000           $-PARAMETER
% IE N                   5000           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SVANBERG';

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
            v_('N') = 10;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('7') = 7;
        v_('8') = 8;
        v_('N-1') = -1+v_('N');
        v_('N-2') = -2+v_('N');
        v_('N-3') = -3+v_('N');
        v_('N-4') = -4+v_('N');
        v_('N-5') = -5+v_('N');
        v_('N-6') = -6+v_('N');
        v_('N-7') = -7+v_('N');
        v_('RN') = v_('N');
        v_('1/N') = 1.0/v_('RN');
        v_('2/N') = 2.0/v_('RN');
        v_('-3/N') = -3.0/v_('RN');
        v_('5/N') = 5.0/v_('RN');
        for I=v_('1'):v_('N')
            v_('RI') = I;
            v_(['B',int2str(I)]) = v_('RI')*v_('5/N');
            v_(['B',int2str(I)]) = 10.0+v_(['B',int2str(I)]);
        end
        for I=v_('1'):v_('2'):v_('N-1')
            v_('RI') = I;
            v_(['A',int2str(I)]) = v_('RI')*v_('2/N');
            v_(['A',int2str(I)]) = 1.0+v_(['A',int2str(I)]);
            v_('I+1') = 1+I;
            v_('RI+1') = v_('I+1');
            v_(['A',int2str(round(v_('I+1')))]) = v_('RI+1')*v_('-3/N');
            v_(['A',int2str(round(v_('I+1')))]) =...
                  5.0+v_(['A',int2str(round(v_('I+1')))]);
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(I)];
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
            pbm.gconst(ig_(['C',int2str(I)])) = v_(['B',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -0.8*ones(pb.n,1);
        pb.xupper = 0.8*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEP',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eEM',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['Q',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEP';
            ielftype(ie) = iet_('eEP');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.8,0.8,[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['P',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEM';
            ielftype(ie) = iet_('eEM');
            vname = ['X',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,-0.8,0.8,[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('2'):v_('N-1')
            v_('I+1') = 1+I;
            ig = ig_(['O',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['A',int2str(I)]);
            ig = ig_(['O',int2str(round(v_('I+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['A',int2str(round(v_('I+1')))]);
        end
        ig = ig_(['C',int2str(round(v_('1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('6')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('6')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('7')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('4')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('6')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('7')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('8')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        for I=v_('5'):v_('2'):v_('N-5')
            v_('I-4') = -4+I;
            v_('I-3') = -3+I;
            v_('I-2') = -2+I;
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            v_('I+2') = 2+I;
            v_('I+3') = 3+I;
            v_('I+4') = 4+I;
            v_('I+5') = 5+I;
            ig = ig_(['C',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I-4')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I-3')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I-2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I-1')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I+1')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I+2')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I+3')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I+4')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['C',int2str(round(v_('I+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I-3')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I-2')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I+2')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I+3')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('I+4')))]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('I+5')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        ig = ig_(['C',int2str(round(v_('N-3')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-7')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-6')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('N-2')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-6')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-5')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-3')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('N-1')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-5')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N')))]);
        pbm.grelw{ig}(posel) = 1.;
        ig = ig_(['C',int2str(round(v_('N')))]);
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('1')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('2')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('3')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('4')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-4')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-3')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N-2')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_(['P',int2str(round(v_('N-1')))]);
        pbm.grelw{ig}(posel) = 1.;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_(['Q',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.;
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN(10)           15.7315
% LO SOLTN(20)           32.4279
% LO SOLTN(30)           49.1425
% LO SOLTN(40)           65.8611
% LO SOLTN(50)           82.5819
% LO SOLTN(60)           99.3039
% LO SOLTN(70)           116.0266
% LO SOLTN(80)           132.7498
% LO SOLTN(90)           149.4734
% LO SOLTN(100)          166.1972
% LO SOLTN(500)          ???
% LO SOLTN(1000)         ???
% LO SOLTN(5000)         ???
% LO SOLTN(10000)        ???
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-V-V';
        pb.x0          = zeros(pb.n,1);
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

    case 'eEP'

        EV_  = varargin{1};
        iel_ = varargin{2};
        D = 1.0/(1.0+EV_(1));
        DSQ = D*D;
        varargout{1} = D;
        if(nargout>1)
            g_(1,1) = -DSQ;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*D*DSQ;
                varargout{3} = H_;
            end
        end

    case 'eEM'

        EV_  = varargin{1};
        iel_ = varargin{2};
        D = 1.0/(1.0-EV_(1));
        DSQ = D*D;
        varargout{1} = D;
        if(nargout>1)
            g_(1,1) = DSQ;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0*D*DSQ;
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

