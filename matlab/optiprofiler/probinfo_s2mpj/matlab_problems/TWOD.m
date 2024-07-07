function varargout = TWOD(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TWOD
%    *********
% 
%    The twod_0 & _00.mod AMPL models from Hans Mittelmann (mittelmann@asu.edu)
%    See: http://plato.asu.edu/ftp/barrier/
% 
%    SIF input: Nick Gould, April 25th 2012
% 
%    classification = 'QLR2-AN-V-V'
% 
%    the x-y discretization 
% 
%       Alternative values for the SIF file parameters:
% IE N                    2             $-PARAMETER
% IE N                   40             $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TWOD';

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
            v_('N') = 2;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   79             $-PARAMETER     twod_000.mod value
% IE N                   99             $-PARAMETER     twod_0.mod value
        v_('M') = v_('N');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('ONE') = 1.0;
        v_('HALF') = 0.5;
        v_('-HALF') = -0.5;
        v_('A') = 0.001;
        v_('UA') = 2.0;
        v_('N1') = -1+v_('N');
        v_('N2') = -2+v_('N');
        v_('M1') = -1+v_('M');
        v_('RN') = v_('N');
        v_('RM') = v_('M');
        v_('DX') = v_('ONE')/v_('RN');
        v_('DY') = v_('ONE')/v_('RM');
        v_('T') = v_('ONE');
        v_('DT') = v_('T')/v_('RM');
        v_('H2') = v_('DX')*v_('DX');
        v_('DXDY') = v_('DX')*v_('DY');
        v_('.5DXDY') = 0.5*v_('DXDY');
        v_('.25DXDY') = 0.25*v_('DXDY');
        v_('.125DXDY') = 0.125*v_('DXDY');
        v_('DTDX') = v_('DT')*v_('DX');
        v_('ADTDX') = v_('A')*v_('DTDX');
        v_('.5ADTDX') = 0.5*v_('ADTDX');
        v_('.25ADTDX') = 0.5*v_('ADTDX');
        v_('1/2DX') = v_('HALF')/v_('DX');
        v_('3/2DX') = 3.0*v_('1/2DX');
        v_('-2/DX') = -4.0*v_('1/2DX');
        v_('3/2DX+1') = 1.0+v_('3/2DX');
        v_('1/2DY') = v_('HALF')/v_('DY');
        v_('3/2DY') = 3.0*v_('1/2DY');
        v_('-2/DY') = -4.0*v_('1/2DY');
        v_('3/2DY+1') = 1.0+v_('3/2DY');
        v_('1/DT') = v_('ONE')/v_('DT');
        v_('-1/DT') = -1.0*v_('1/DT');
        v_('-.1/2H2') = v_('-HALF')/v_('H2');
        v_('2/H2') = -4.0*v_('-.1/2H2');
        v_('1/DT+2/H2') = v_('1/DT')+v_('2/H2');
        v_('-1/DT+2/H2') = v_('-1/DT')+v_('2/H2');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                for K=v_('0'):v_('M')
                    [iv,ix_] =...
                          s2mpjlib('ii',['Y',int2str(K),',',int2str(I),',',int2str(J)],ix_);
                    pb.xnames{iv} = ['Y',int2str(K),',',int2str(I),',',int2str(J)];
                end
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('0'):v_('N1')
                [iv,ix_] = s2mpjlib('ii',['U',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['U',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N1')
            v_('I+') = 1+I;
            v_('I-') = -1+I;
            for J=v_('1'):v_('N1')
                v_('J+') = 1+J;
                v_('J-') = -1+J;
                for K=v_('0'):v_('M1')
                    v_('K+') = 1+K;
                    [ig,ig_] =...
                          s2mpjlib('ii',['P',int2str(K),',',int2str(I),',',int2str(J)],ig_);
                    gtype{ig}  = '==';
                    cnames{ig} = ['P',int2str(K),',',int2str(I),',',int2str(J)];
                    iv = ix_(['Y',int2str(round(v_('K+'))),',',int2str(I),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('1/DT+2/H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('1/DT+2/H2');
                    end
                    iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-1/DT+2/H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-1/DT+2/H2');
                    end
                    iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('J-')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('J+')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv = ix_(['Y',int2str(K),',',int2str(round(v_('I-'))),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv = ix_(['Y',int2str(K),',',int2str(round(v_('I+'))),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv =...
                          ix_(['Y',int2str(round(v_('K+'))),',',int2str(round(v_('I-'))),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv =...
                          ix_(['Y',int2str(round(v_('K+'))),',',int2str(round(v_('I+'))),',',int2str(J)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv =...
                          ix_(['Y',int2str(round(v_('K+'))),',',int2str(I),',',int2str(round(v_('J-')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                    iv =...
                          ix_(['Y',int2str(round(v_('K+'))),',',int2str(I),',',int2str(round(v_('J+')))]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-.1/2H2')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-.1/2H2');
                    end
                end
            end
        end
        for I=v_('1'):v_('N1')
            for K=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['B1',int2str(K),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['B1',int2str(K),',',int2str(I)];
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('N2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/2DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/2DY');
                end
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('N1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/DY');
                end
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('N')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('3/2DY+1')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('3/2DY+1');
                end
                iv = ix_(['U',int2str(K),',',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
                [ig,ig_] = s2mpjlib('ii',['B2',int2str(K),',',int2str(I)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['B2',int2str(K),',',int2str(I)];
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('2')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/2DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/2DY');
                end
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('1')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/DY')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/DY');
                end
                iv = ix_(['Y',int2str(K),',',int2str(I),',',int2str(round(v_('0')))]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('3/2DY+1')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('3/2DY+1');
                end
            end
        end
        for J=v_('1'):v_('N1')
            for K=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['B3',int2str(K),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['B3',int2str(K),',',int2str(J)];
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('2'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/2DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/2DX');
                end
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/DX');
                end
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('0'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('3/2DX+1')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('3/2DX+1');
                end
                [ig,ig_] = s2mpjlib('ii',['B4',int2str(K),',',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['B4',int2str(K),',',int2str(J)];
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('N2'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('1/2DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('1/2DX');
                end
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('N1'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('-2/DX')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('-2/DX');
                end
                iv = ix_(['Y',int2str(K),',',int2str(round(v_('N'))),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_('3/2DX+1')+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_('3/2DX+1');
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
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                pb.xlower(ix_(['Y',int2str(round(v_('0'))),',',int2str(I),',',int2str(J)]),1) = 0.0;
                pb.xupper(ix_(['Y',int2str(round(v_('0'))),',',int2str(I),',',int2str(J)]),1) = 0.0;
            end
        end
        for I=v_('0'):v_('N')
            for J=v_('0'):v_('N')
                for K=v_('1'):v_('M')
                    pb.xlower(ix_(['Y',int2str(K),',',int2str(I),',',int2str(J)]),1) = 0.0;
                    pb.xupper(ix_(['Y',int2str(K),',',int2str(I),',',int2str(J)])) = 0.8;
                end
            end
        end
        for I=v_('1'):v_('M')
            for J=v_('0'):v_('N1')
                pb.xlower(ix_(['U',int2str(I),',',int2str(J)]),1) = 0.0;
                pb.xupper(ix_(['U',int2str(I),',',int2str(J)])) = v_('UA');
            end
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        pb.y0 = 0.0*ones(pb.m,1);
        for I=v_('1'):v_('M')
            for J=v_('0'):v_('N1')
                if(isKey(ix_,['U',int2str(I),',',int2str(J)]))
                    pb.x0(ix_(['U',int2str(I),',',int2str(J)]),1) = v_('UA');
                else
                    pb.y0(find(pbm.congrps==ig_(['U',int2str(I),',',int2str(J)])),1) =...
                          v_('UA');
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQ',iet_);
        elftv{it}{1} = 'U';
        [it,iet_] = s2mpjlib( 'ii', 'eSQD',iet_);
        elftv{it}{1} = 'Y';
        elftp{it}{1} = 'YP';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('0'):v_('N')
            v_('RI') = I;
            v_('.5DXDYI') = v_('.5DXDY')*v_('RI');
            for J=v_('0'):v_('N')
                v_('RJ') = J;
                v_('.5DXDYIJ') = v_('.5DXDYI')*v_('RJ');
                v_('YP') = 0.25+v_('.5DXDYIJ');
                ename = ['E',int2str(round(v_('M'))),',',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQD';
                ielftype(ie) = iet_('eSQD');
                ename = ['E',int2str(round(v_('M'))),',',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                vname = ['Y',int2str(round(v_('M'))),',',int2str(I),',',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                ename = ['E',int2str(round(v_('M'))),',',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                [~,posep] = ismember('YP',elftp{ielftype(ie)});
                pbm.elpar{ie}(posep) = v_('YP');
            end
        end
        for K=v_('1'):v_('M')
            for I=v_('1'):v_('N1')
                ename = ['E',int2str(K),',',int2str(I)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'eSQ';
                ielftype(ie) = iet_('eSQ');
                vname = ['U',int2str(K),',',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
                posev = find(strcmp('U',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        ig = ig_('OBJ');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('0'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('.125DXDY');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('0'))),',',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('.125DXDY');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('N'))),',',int2str(round(v_('0')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('.125DXDY');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) =...
              ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('.125DXDY');
        for J=v_('1'):v_('N1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) =...
                  ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('0'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('.25DXDY');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) =...
                  ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('N'))),',',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('.25DXDY');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) =...
                  ie_(['E',int2str(round(v_('M'))),',',int2str(J),',',int2str(round(v_('0')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('.25DXDY');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) =...
                  ie_(['E',int2str(round(v_('M'))),',',int2str(round(v_('N'))),',',int2str(round(v_('N')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('.25DXDY');
        end
        for I=v_('1'):v_('N1')
            for J=v_('1'):v_('N1')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['E',int2str(round(v_('M'))),',',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('.5DXDY');
            end
        end
        for K=v_('1'):v_('M1')
            for I=v_('1'):v_('N1')
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(K),',',int2str(I)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('.5ADTDX');
            end
        end
        for I=v_('1'):v_('N1')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('M'))),',',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('.25ADTDX');
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'QLR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end

    case 'eSQD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = (EV_(1)-pbm.elpar{iel_}(1))^2;
        if(nargout>1)
            g_(1,1) = 2.0*(EV_(1)-pbm.elpar{iel_}(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
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

