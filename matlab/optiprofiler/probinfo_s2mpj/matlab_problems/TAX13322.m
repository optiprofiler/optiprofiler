function varargout = TAX13322(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : TAX13322
%    --------
% 
%    An optimal income tax model with multidimensional taxpayer types,
%    due to Judd, Ma, Saunders & Su
% 
%    Source:
%    Kenneth L. Judd, Ma,  Michael A. Saunders and Che-Lin Su
%    "Optimal Income Taxation with Multidimensional Taxpayer Types"
%    Working Paper, Hoover Institute, Stanford University, 2017
% 
%    SIF input: Nick Gould, July 2018 based on the AMPL model pTAX5Dncl
% 
%    "If ever there was an example that exhibited the stupidity of SIF,
%     this is it. NIMG"
% 
%    classification = 'OOR2-MN-72-1261'
% 
%    parameters
% 
%       Alternative values for the SIF file parameters:
% IE NA                  1              $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'TAX13322';

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
            v_('NA') = 1;  %  SIF file default value
        else
            v_('NA') = varargin{1};
        end
% IE NB                  3              $-PARAMETER
        if(nargs<2)
            v_('NB') = 3;  %  SIF file default value
        else
            v_('NB') = varargin{2};
        end
% IE NC                  3              $-PARAMETER
        if(nargs<3)
            v_('NC') = 3;  %  SIF file default value
        else
            v_('NC') = varargin{3};
        end
% IE ND                  2              $-PARAMETER
        if(nargs<4)
            v_('ND') = 2;  %  SIF file default value
        else
            v_('ND') = varargin{4};
        end
% IE NE                  2              $-PARAMETER
        if(nargs<5)
            v_('NE') = 2;  %  SIF file default value
        else
            v_('NE') = varargin{5};
        end
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('6') = 6;
        v_('ONE') = 1.0e0;
        v_('TWO') = 2.0e0;
        v_('THREE') = 3.0e0;
        v_('NBD') = v_('NB')*v_('ND');
        v_('NCE') = v_('NC')*v_('NE');
        v_('NP') = v_('NBD')*v_('NCE');
        v_('NP') = v_('NP')*v_('NA');
        v_('NPM1') = -1+v_('NP');
        v_('M') = v_('NP')*v_('NPM1');
        v_('OMEGA1') = v_('ONE')/v_('TWO');
        v_('OMEGA2') = v_('TWO')/v_('THREE');
        v_('THETA1') = v_('ONE')/v_('THREE');
        v_('THETA2') = v_('ONE')/v_('TWO');
        v_('THETA3') = v_('TWO')/v_('THREE');
        v_('PSI1') = 1.0e0;
        v_('PSI2') = 1.5e0;
        v_('W1') = 2.0e0;
        v_('W2') = 2.5e0;
        v_('W3') = 3.0e0;
        v_('W4') = 3.5e0;
        v_('W5') = 4.0e0;
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    v_(['LAM',int2str(I),',',int2str(P),',',int2str(Q)]) = 1.0e+0;
                end
            end
        end
        v_('Q') = v_('1');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA1');
        v_('Q') = v_('2');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA2');
        v_('Q') = v_('3');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA1');
        v_('Q') = v_('4');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA2');
        v_('Q') = v_('5');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA1');
        v_('Q') = v_('6');
        v_(['RA',int2str(round(v_('Q')))]) = v_('ONE')/v_('OMEGA2');
        for I=v_('1'):v_('NA')
            v_('LOGW') = log(v_(['W',int2str(I)]));
            v_('P') = v_('1');
            v_('-THETA') = -1.0*v_('THETA1');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI1');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
            v_('P') = v_('2');
            v_('-THETA') = -1.0*v_('THETA1');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI2');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
            v_('P') = v_('3');
            v_('-THETA') = -1.0*v_('THETA2');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI1');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
            v_('P') = v_('4');
            v_('-THETA') = -1.0*v_('THETA2');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI2');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
            v_('P') = v_('5');
            v_('-THETA') = -1.0*v_('THETA3');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI1');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
            v_('P') = v_('6');
            v_('-THETA') = -1.0*v_('THETA3');
            v_('RB') = v_('LOGW')*v_('-THETA');
            v_('RE') = exp(v_('RB'));
            v_('RB') = v_('RE')*v_('PSI2');
            v_(['RB',int2str(I),',',int2str(round(v_('P')))]) = v_('RB')/v_('-THETA');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    [iv,ix_] =...
                          s2mpjlib('ii',['C',int2str(I),',',int2str(P),',',int2str(Q)],ix_);
                    pb.xnames{iv} = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_] =...
                          s2mpjlib('ii',['Y',int2str(I),',',int2str(P),',',int2str(Q)],ix_);
                    pb.xnames{iv} = ['Y',int2str(I),',',int2str(P),',',int2str(Q)];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for L=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['I',int2str(L)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['I',int2str(L)];
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    v_('LAMBDA') = v_(['LAM',int2str(I),',',int2str(P),',',int2str(Q)]);
                    v_('-LAMBDA') = -1.0e0*v_('LAMBDA');
                    [ig,ig_] = s2mpjlib('ii','T',ig_);
                    gtype{ig}  = '>=';
                    cnames{ig} = 'T';
                    iv = ix_(['Y',int2str(I),',',int2str(P),',',int2str(Q)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('LAMBDA')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('LAMBDA');
                    end
                    iv = ix_(['C',int2str(I),',',int2str(P),',',int2str(Q)]);
                    if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                        pbm.A(ig,iv) = v_('-LAMBDA')+pbm.A(ig,iv);
                    else
                        pbm.A(ig,iv) = v_('-LAMBDA');
                    end
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
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    pb.xlower(ix_(['C',int2str(I),',',int2str(P),',',int2str(Q)]),1) = 0.1e0;
                    pb.xlower(ix_(['Y',int2str(I),',',int2str(P),',',int2str(Q)]),1) = 0.1e0;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.1e0*ones(pb.n,1);
        pb.y0 = 0.1e0*ones(pb.m,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eA1',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eA2',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eA3',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eA4',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eA5',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eA6',iet_);
        elftv{it}{1} = 'C';
        [it,iet_] = s2mpjlib( 'ii', 'eB1',iet_);
        elftv{it}{1} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eB2',iet_);
        elftv{it}{1} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eB3',iet_);
        elftv{it}{1} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A1-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA1';
                    ielftype(ie) = iet_('eA1');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A2-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA2';
                    ielftype(ie) = iet_('eA2');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A3-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA3';
                    ielftype(ie) = iet_('eA3');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A4-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA4';
                    ielftype(ie) = iet_('eA4');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A5-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA5';
                    ielftype(ie) = iet_('eA5');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['A6-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eA6';
                    ielftype(ie) = iet_('eA6');
                    vname = ['C',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('C',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['B1-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eB1';
                    ielftype(ie) = iet_('eB1');
                    vname = ['Y',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('Y',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['B2-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eB2';
                    ielftype(ie) = iet_('eB2');
                    vname = ['Y',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('Y',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                for Q=v_('1'):v_('NCE')
                    ename = ['B3-',int2str(I),',',int2str(P),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eB3';
                    ielftype(ie) = iet_('eB3');
                    vname = ['Y',int2str(I),',',int2str(P),',',int2str(Q)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.1e0);
                    posev = find(strcmp('Y',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('NA')
            for P=v_('1'):v_('NBD')
                v_('Q') = v_('1');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('Q') = v_('2');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('Q') = v_('3');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('Q') = v_('4');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('Q') = v_('5');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('Q') = v_('6');
                v_('R') =...
                      v_(['RA',int2str(round(v_('Q')))])*v_(['LAM',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(P),',',int2str(round(v_('Q')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
            end
        end
        for I=v_('1'):v_('NA')
            v_('LOGW') = log(v_(['W',int2str(I)]));
            for Q=v_('1'):v_('NCE')
                v_('P') = v_('1');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                ig = ig_('OBJ');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('P') = v_('2');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('P') = v_('3');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('P') = v_('4');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('P') = v_('5');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
                v_('P') = v_('6');
                v_('R') =...
                      v_(['RB',int2str(I),',',int2str(round(v_('P')))])*v_(['LAM',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('R');
            end
        end
        v_('L') = v_('0');
        for I=v_('1'):v_('NA')
            v_('I+1') = 1+I;
            v_('I-1') = -1+I;
            v_('LOGW') = log(v_(['W',int2str(I)]));
            v_('P') = v_('1');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('P') = v_('2');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('P') = v_('3');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('P') = v_('4');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('P') = v_('5');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('P') = v_('6');
            v_('P+1') = 1+v_('P');
            v_('P-1') = -1+v_('P');
            v_('RB') = v_(['RB',int2str(I),',',int2str(round(v_('P')))]);
            v_('-RB') = -1.0*v_('RB');
            v_('RA') = v_(['RA',int2str(round(v_('1')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A1-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('2'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A1-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('1')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('2')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A2-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('1')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('3'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A2-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('2')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('3')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('2')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('4'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('3')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('4')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A4-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('3')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('5'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A4-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('4')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('5')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A5-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('4')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            for T=v_('6'):v_('NCE')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A5-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('5')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
            v_('RA') = v_(['RA',int2str(round(v_('6')))]);
            v_('-RA') = -1.0*v_('RA');
            for S=v_('1'):v_('P-1')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for S=v_('P+1'):v_('NBD')
                for R=v_('1'):v_('NA')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('1')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('2')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('3')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('4')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('5')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                    v_('L') = v_('L')+v_('1');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(S),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('1'):v_('I-1')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for R=v_('I+1'):v_('NA')
                for T=v_('1'):v_('NCE')
                    v_('L') = v_('L')+v_('1');
                    ig = ig_(['I',int2str(round(v_('L')))]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['A6-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RA');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('RB');
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) =...
                          ie_(['B3-',int2str(R),',',int2str(round(v_('P'))),',',int2str(T)]);
                    nlc = union(nlc,ig);
                    pbm.grelw{ig}(posel) = v_('-RB');
                end
            end
            for T=v_('1'):v_('5')
                v_('L') = v_('L')+v_('1');
                ig = ig_(['I',int2str(round(v_('L')))]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['A6-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RA');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(round(v_('6')))]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('RB');
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) =...
                      ie_(['B3-',int2str(I),',',int2str(round(v_('P'))),',',int2str(T)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-RB');
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'OOR2-MN-72-1261';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 0.1e0;
        varargout{1} = pbm;

    case 'eA1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 0.0e0;
        OMEGA = 1.0e0/2.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eA2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 0.0e0;
        OMEGA = 2.0e0/3.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eA3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 1.0e0;
        OMEGA = 1.0e0/2.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eA4'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 1.0e0;
        OMEGA = 2.0e0/3.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eA5'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 1.5e0;
        OMEGA = 1.0e0/2.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eA6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ALPHA = 1.5e0;
        OMEGA = 2.0e0/3.0e0;
        OM1 = OMEGA-1.0e0;
        OM2 = OMEGA-2.0e0;
        CMA = EV_(1)-ALPHA;
        BIG = CMA>=pbm.efpar(1);
        if(BIG)
            F = CMA^OMEGA;
        end
        if(BIG)
            G = OMEGA*CMA^OM1;
        end
        if(BIG)
            H = OMEGA*OM1*CMA^OM2;
        end
        if(~BIG)
            C1 = OMEGA*(2.0e0-OMEGA)*pbm.efpar(1)^OM1;
        end
        if(~BIG)
            C2 = OMEGA*OM1*pbm.efpar(1)^OM2;
        end
        if(~BIG)
            F = (1.0e0-1.5e0*OMEGA+0.5e0*OMEGA*OMEGA)*pbm.efpar(1)^OMEGA+...
             C1*CMA+0.5e0*C2*CMA*CMA;
        end
        if(~BIG)
            G = C1+C2*CMA;
        end
        if(~BIG)
            H = C2;
        end
        varargout{1} = F;
        if(nargout>1)
            g_(1,1) = G;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = H;
                varargout{3} = H_;
            end
        end

    case 'eB1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        THETA = 1.0e0/3.0e0;
        varargout{1} = EV_(1)^THETA;
        if(nargout>1)
            g_(1,1) = THETA*EV_(1)^(THETA-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = THETA*(THETA-1.0e0)*EV_(1)^(THETA-2.0e0);
                varargout{3} = H_;
            end
        end

    case 'eB2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        THETA = 1.0e0/2.0e0;
        varargout{1} = EV_(1)^THETA;
        if(nargout>1)
            g_(1,1) = THETA*EV_(1)^(THETA-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = THETA*(THETA-1.0e0)*EV_(1)^(THETA-2.0e0);
                varargout{3} = H_;
            end
        end

    case 'eB3'

        EV_  = varargin{1};
        iel_ = varargin{2};
        THETA = 2.0e0/3.0e0;
        varargout{1} = EV_(1)^THETA;
        if(nargout>1)
            g_(1,1) = THETA*EV_(1)^(THETA-1.0e0);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = THETA*(THETA-1.0e0)*EV_(1)^(THETA-2.0e0);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [1,0];
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

