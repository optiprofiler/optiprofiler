function varargout = FERRISDC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : FERRISDC
%    *********
% 
%    A QP suggested by Michael Ferris
%    classification = ''
%    SIF input: Nick Gould, November 2001.
% 
%    classification = 'QLR2-AN-V-V'
% 
%       Alternative values for the SIF file parameters:
% IE n                   4              $-PARAMETER
% IE n                   100            $-PARAMETER
% IE n                   200            $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'FERRISDC';

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
            v_('n') = 4;  %  SIF file default value
        else
            v_('n') = varargin{1};
        end
% IE n                   300            $-PARAMETER
% IE k                   3              $-PARAMETER
% IE k                   10             $-PARAMETER
        if(nargs<2)
            v_('k') = 3;  %  SIF file default value
        else
            v_('k') = varargin{2};
        end
% IE k                   20             $-PARAMETER
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('12') = 12.0;
        v_('24') = 24.0;
        v_('240') = 240.0;
        v_('k-1') = -1+v_('k');
        v_('k') = v_('k');
        v_('k-1') = v_('k-1');
        v_('-1/k-1') = -1.0/v_('k-1');
        v_('-1/k') = -1.0/v_('k');
        v_('n') = v_('n');
        v_('1') = 1.0;
        v_('2') = 2.0;
        v_('1/12') = v_('1')/v_('12');
        v_('1/24') = v_('1')/v_('24');
        v_('1/240') = v_('1')/v_('240');
        v_('7/240') = 7.0*v_('1/240');
        v_('2**2') = v_('2')*v_('2');
        v_('2**4') = v_('2**2')*v_('2**2');
        v_('2**8') = v_('2**4')*v_('2**4');
        v_('2**10') = v_('2**8')*v_('2**2');
        v_('2**16') = v_('2**8')*v_('2**8');
        v_('2**26') = v_('2**16')*v_('2**10');
        v_('2**-26') = v_('1')/v_('2**26');
        v_('nlambda') = v_('n')*v_('2**-26');
        v_('-1/k-1*nl') = v_('nlambda')*v_('-1/k-1');
        v_('ix') = 1;
        v_('ax') = 16807;
        v_('b15') = 32768;
        v_('b16') = 65536;
        v_('pp') = 2147483647;
        v_('pp') = v_('pp');
        for j=v_('1'):v_('n')
            v_('xhi') = fix(v_('ix')/v_('b16'));
            v_('xalo') = v_('xhi')*v_('b16');
            v_('xalo') = v_('ix')-v_('xalo');
            v_('xalo') = v_('xalo')*v_('ax');
            v_('leftlo') = fix(v_('xalo')/v_('b16'));
            v_('fhi') = v_('xhi')*v_('ax');
            v_('fhi') = v_('fhi')+v_('leftlo');
            v_('kk') = fix(v_('fhi')/v_('b15'));
            v_('dum') = v_('leftlo')*v_('b16');
            v_('dum') = v_('xalo')-v_('dum');
            v_('ix') = v_('dum')-v_('pp');
            v_('dum') = v_('kk')*v_('b15');
            v_('dum') = v_('fhi')-v_('dum');
            v_('dum') = v_('dum')*v_('b16');
            v_('ix') = v_('ix')+v_('dum');
            v_('ix') = v_('ix')+v_('kk');
            v_('a') = v_('ix');
            v_('a') = -1.0*v_('a');
            v_('b') = 0.0;
            v_('absa') = abs(v_('a'));
            v_('absb') = abs(v_('b'));
            v_('absa+b') = v_('absa')+v_('absb');
            v_('absa+b+2') = 2.0+v_('absa+b');
            v_('a') = v_('a')+v_('absa+b+2');
            v_('b') = v_('b')+v_('absa+b+2');
            v_('a/b') = v_('a')/v_('b');
            v_('b/a') = v_('b')/v_('a');
            v_('a/b') = fix(v_('a/b'));
            v_('b/a') = fix(v_('b/a'));
            v_('a/b') = v_('a/b');
            v_('b/a') = v_('b/a');
            v_('sum') = v_('a/b')+v_('b/a');
            v_('a') = v_('a')*v_('a/b');
            v_('b') = v_('b')*v_('b/a');
            v_('maxa,b') = v_('a')+v_('b');
            v_('maxa,b') = v_('maxa,b')/v_('sum');
            v_('c') = v_('absa+b+2')-v_('maxa,b');
            v_('a') = v_('absa+b+2')-v_('a');
            v_('absc') = abs(v_('c'));
            v_('absc+1') = 1.0+v_('absc');
            v_('absc+2') = 2.0+v_('absc');
            v_('f') = v_('absc+2')/v_('absc+1');
            v_('f') = fix(v_('f'));
            v_('g') = 2-v_('f');
            for l=v_('1'):v_('g')
                v_('ix') = v_('ix')+v_('pp');
            end
            v_('randp') = v_('ix');
            v_(['X',int2str(j)]) = v_('randp')/v_('pp');
        end
        for j=v_('1'):v_('n')
            v_('xhi') = fix(v_('ix')/v_('b16'));
            v_('xalo') = v_('xhi')*v_('b16');
            v_('xalo') = v_('ix')-v_('xalo');
            v_('xalo') = v_('xalo')*v_('ax');
            v_('leftlo') = fix(v_('xalo')/v_('b16'));
            v_('fhi') = v_('xhi')*v_('ax');
            v_('fhi') = v_('fhi')+v_('leftlo');
            v_('kk') = fix(v_('fhi')/v_('b15'));
            v_('dum') = v_('leftlo')*v_('b16');
            v_('dum') = v_('xalo')-v_('dum');
            v_('ix') = v_('dum')-v_('pp');
            v_('dum') = v_('kk')*v_('b15');
            v_('dum') = v_('fhi')-v_('dum');
            v_('dum') = v_('dum')*v_('b16');
            v_('ix') = v_('ix')+v_('dum');
            v_('ix') = v_('ix')+v_('kk');
            v_('a') = v_('ix');
            v_('a') = -1.0*v_('a');
            v_('b') = 0.0;
            v_('absa') = abs(v_('a'));
            v_('absb') = abs(v_('b'));
            v_('absa+b') = v_('absa')+v_('absb');
            v_('absa+b+2') = 2.0+v_('absa+b');
            v_('a') = v_('a')+v_('absa+b+2');
            v_('b') = v_('b')+v_('absa+b+2');
            v_('a/b') = v_('a')/v_('b');
            v_('b/a') = v_('b')/v_('a');
            v_('a/b') = fix(v_('a/b'));
            v_('b/a') = fix(v_('b/a'));
            v_('a/b') = v_('a/b');
            v_('b/a') = v_('b/a');
            v_('sum') = v_('a/b')+v_('b/a');
            v_('a') = v_('a')*v_('a/b');
            v_('b') = v_('b')*v_('b/a');
            v_('maxa,b') = v_('a')+v_('b');
            v_('maxa,b') = v_('maxa,b')/v_('sum');
            v_('c') = v_('absa+b+2')-v_('maxa,b');
            v_('a') = v_('absa+b+2')-v_('a');
            v_('absc') = abs(v_('c'));
            v_('absc+1') = 1.0+v_('absc');
            v_('absc+2') = 2.0+v_('absc');
            v_('f') = v_('absc+2')/v_('absc+1');
            v_('f') = fix(v_('f'));
            v_('g') = 2-v_('f');
            for l=v_('1'):v_('g')
                v_('ix') = v_('ix')+v_('pp');
            end
            v_('randp') = v_('ix');
            v_(['R',int2str(j)]) = v_('randp')/v_('pp');
        end
        for j=v_('1'):v_('n')
            v_('arg') = -3.0*v_(['X',int2str(j)]);
            v_('arg') = exp(v_('arg'));
            v_(['P',int2str(j),',',int2str(round(v_('1')))]) = 0.97*v_('arg');
            v_('arg') = -1.2+v_(['X',int2str(j)]);
            v_('arg') = v_('arg')*v_('arg');
            v_('arg') = -2.5*v_('arg');
            v_(['P',int2str(j),',',int2str(round(v_('3')))]) = exp(v_('arg'));
            v_('arg') = v_('1')-v_(['P',int2str(j),',',int2str(round(v_('1')))]);
            v_(['P',int2str(j),',',int2str(round(v_('2')))]) =...
                  v_('arg')-v_(['P',int2str(j),',',int2str(round(v_('3')))]);
        end
        for j=v_('1'):v_('n')
            v_('a') =...
                  v_(['P',int2str(j),',',int2str(round(v_('1')))])-v_(['R',int2str(j)]);
            v_('a') = -1.0*v_('a');
            v_('b') = 0.0;
            v_('absa') = abs(v_('a'));
            v_('absb') = abs(v_('b'));
            v_('absa+b') = v_('absa')+v_('absb');
            v_('absa+b+2') = 2.0+v_('absa+b');
            v_('a') = v_('a')+v_('absa+b+2');
            v_('b') = v_('b')+v_('absa+b+2');
            v_('a/b') = v_('a')/v_('b');
            v_('b/a') = v_('b')/v_('a');
            v_('a/b') = fix(v_('a/b'));
            v_('b/a') = fix(v_('b/a'));
            v_('a/b') = v_('a/b');
            v_('b/a') = v_('b/a');
            v_('sum') = v_('a/b')+v_('b/a');
            v_('a') = v_('a')*v_('a/b');
            v_('b') = v_('b')*v_('b/a');
            v_('maxa,b') = v_('a')+v_('b');
            v_('maxa,b') = v_('maxa,b')/v_('sum');
            v_('c') = v_('absa+b+2')-v_('maxa,b');
            v_('a') = v_('absa+b+2')-v_('a');
            v_('absc') = abs(v_('c'));
            v_('absc+1') = 1.0+v_('absc');
            v_('absc+2') = 2.0+v_('absc');
            v_('f') = v_('absc+2')/v_('absc+1');
            v_('f') = fix(v_('f'));
            v_('g') = 2-v_('f');
            for l1=v_('g'):v_('0')
                v_(['y',int2str(j)]) = 1.0;
            end
            for l1=v_('1'):v_('g')
                v_('a') = v_('1')-v_(['P',int2str(j),',',int2str(round(v_('3')))]);
                v_('a') = v_('a')-v_(['R',int2str(j)]);
                v_('a') = -1.0*v_('a');
                v_('b') = 0.0;
                v_('absa') = abs(v_('a'));
                v_('absb') = abs(v_('b'));
                v_('absa+b') = v_('absa')+v_('absb');
                v_('absa+b+2') = 2.0+v_('absa+b');
                v_('a') = v_('a')+v_('absa+b+2');
                v_('b') = v_('b')+v_('absa+b+2');
                v_('a/b') = v_('a')/v_('b');
                v_('b/a') = v_('b')/v_('a');
                v_('a/b') = fix(v_('a/b'));
                v_('b/a') = fix(v_('b/a'));
                v_('a/b') = v_('a/b');
                v_('b/a') = v_('b/a');
                v_('sum') = v_('a/b')+v_('b/a');
                v_('a') = v_('a')*v_('a/b');
                v_('b') = v_('b')*v_('b/a');
                v_('maxa,b') = v_('a')+v_('b');
                v_('maxa,b') = v_('maxa,b')/v_('sum');
                v_('c') = v_('absa+b+2')-v_('maxa,b');
                v_('a') = v_('absa+b+2')-v_('a');
                v_('absc') = abs(v_('c'));
                v_('absc+1') = 1.0+v_('absc');
                v_('absc+2') = 2.0+v_('absc');
                v_('f') = v_('absc+2')/v_('absc+1');
                v_('f') = fix(v_('f'));
                v_('g') = 2-v_('f');
                for l2=v_('g'):v_('0')
                    v_(['y',int2str(j)]) = 2.0;
                end
                for l2=v_('1'):v_('g')
                    v_(['y',int2str(j)]) = 3.0;
                end
            end
        end
        for j=v_('1'):v_('n')
            v_('yj') = v_(['y',int2str(j)]);
            v_('yj') = fix(v_('yj'));
            for i=v_('1'):v_('k')
                v_('c') = v_('yj')-i;
                v_('c') = v_('c');
                v_('absc') = abs(v_('c'));
                v_('absc+1') = 1.0+v_('absc');
                v_('absc+2') = 2.0+v_('absc');
                v_('f') = v_('absc+2')/v_('absc+1');
                v_('f') = fix(v_('f'));
                v_('g') = 2-v_('f');
                for l=v_('g'):v_('0')
                    v_(['Y',int2str(i),',',int2str(j)]) = v_('nlambda');
                end
                for l=v_('1'):v_('g')
                    v_(['Y',int2str(i),',',int2str(j)]) = v_('-1/k-1*nl');
                end
            end
        end
        for i=v_('1'):v_('n')
            v_('di') = -0.5+v_(['X',int2str(i)]);
            v_('di2') = v_('di')*v_('di');
            v_('di2') = v_('di2')-v_('1/12');
            for j=i:v_('n')
                v_('Xi-Xj') = v_(['X',int2str(i)])-v_(['X',int2str(j)]);
                v_('bij') = abs(v_('Xi-Xj'));
                v_('dj') = -0.5+v_(['X',int2str(j)]);
                v_('dj2') = v_('dj')*v_('dj');
                v_('dj2') = v_('dj2')-v_('1/12');
                v_('c') = -0.5+v_('bij');
                v_('c2') = v_('c')*v_('c');
                v_('c4') = v_('c2')*v_('c2');
                v_('c2') = -0.5*v_('c2');
                v_('arg') = v_('7/240')+v_('c2');
                v_('arg') = v_('arg')+v_('c4');
                v_('arg') = v_('arg')*v_('1/24');
                v_('dij') = v_('di')*v_('dj');
                v_('dij2') = v_('di2')*v_('dj2');
                v_('dij2') = 0.25*v_('dij2');
                v_('arg') = v_('dij2')-v_('arg');
                v_(['K',int2str(i),',',int2str(j)]) = v_('dij')+v_('arg');
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for i=v_('1'):v_('k')
            for j=v_('1'):v_('n')
                [iv,ix_] = s2mpjlib('ii',['A',int2str(i),',',int2str(j)],ix_);
                pb.xnames{iv} = ['A',int2str(i),',',int2str(j)];
            end
        end
        for i=v_('1'):v_('n')
            [iv,ix_] = s2mpjlib('ii',['W',int2str(i)],ix_);
            pb.xnames{iv} = ['W',int2str(i)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for j=v_('1'):v_('n')
            for i=v_('1'):v_('k')
                [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
                gtype{ig} = '<>';
                iv = ix_(['A',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['Y',int2str(i),',',int2str(j)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['Y',int2str(i),',',int2str(j)]);
                end
            end
        end
        for i=v_('1'):v_('k')
            for j=v_('1'):v_('n')
                [ig,ig_] = s2mpjlib('ii',['C',int2str(i)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['C',int2str(i)];
                iv = ix_(['A',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['W',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-1/k')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-1/k');
                end
            end
        end
        for j=v_('1'):v_('n')
            [ig,ig_] = s2mpjlib('ii',['A',int2str(j)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['A',int2str(j)];
            iv = ix_(['W',int2str(j)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
            for i=v_('1'):v_('k')
                [ig,ig_] = s2mpjlib('ii',['A',int2str(j)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['A',int2str(j)];
                iv = ix_(['A',int2str(i),',',int2str(j)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
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
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = 0.0*ones(pb.n,1);
        pb.xupper = 1.0*ones(pb.n,1);
        for i=v_('1'):v_('n')
            pb.xlower(ix_(['W',int2str(i)])) = -Inf;
            pb.xupper(ix_(['W',int2str(i)]),1) = +Inf;
        end
        for j=v_('1'):v_('n')
            v_('yj') = v_(['y',int2str(j)]);
            v_('yj') = fix(v_('yj'));
            for i=v_('1'):v_('k')
                v_('c') = v_('yj')-i;
                v_('c') = v_('c');
                v_('absc') = abs(v_('c'));
                v_('absc+1') = 1.0+v_('absc');
                v_('absc+2') = 2.0+v_('absc');
                v_('f') = v_('absc+2')/v_('absc+1');
                v_('f') = fix(v_('f'));
                v_('g') = 2-v_('f');
                for l=v_('g'):v_('0')
                    pb.xlower(ix_(['A',int2str(i),',',int2str(j)]),1) = 0.0;
                    pb.xupper(ix_(['A',int2str(i),',',int2str(j)]),1) = 0.0;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% QUADRATIC %%%%%%%%%%%%%%%%%%%
        pbm.H = sparse( pb.n, pb.n );
        for i=v_('1'):v_('k')
            for l=v_('1'):v_('n')
                for j=v_('1'):l
                    ix1 = ix_(['A',int2str(i),',',int2str(j)]);
                    ix2 = ix_(['A',int2str(i),',',int2str(l)]);
                    pbm.H(ix1,ix2) = v_(['K',int2str(j),',',int2str(l)])+pbm.H(ix1,ix2);
                    pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
                end
            end
        end
        for l=v_('1'):v_('n')
            for j=v_('1'):l
                v_('coef') = v_('-1/k')*v_(['K',int2str(j),',',int2str(l)]);
                ix1 = ix_(['W',int2str(j)]);
                ix2 = ix_(['W',int2str(l)]);
                pbm.H(ix1,ix2) = v_('coef')+pbm.H(ix1,ix2);
                pbm.H(ix2,ix1) = pbm.H(ix1,ix2);
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% XL SOLUTION            -1.131846D+2   $ nlambda = 1.5625
% XL SOLUTION            -8.032841E-5   $ nlambda = 1.4901E-06
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.lincons   = [1:length(pbm.congrps)];
        pb.pbclass = 'QLR2-AN-V-V';
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

