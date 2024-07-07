function varargout = NELSON(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : NELSON
%    *********
% 
%    NIST Data fitting problem NELSON given as an inconsistent set of
%    nonlinear equations.
% 
%    Fit: log[y] = b1 - b2*x1 * exp[-b3*x2] + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    Reference: Nelson, W. (1981).  
%      Analysis of Performance-Degradation Data.  
%      IEEE Transactions on Reliability. Vol. 2, R-30, No. 2, pp. 149-155.
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%    classification = 'NOR2-MN-3-128'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'NELSON';

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
        v_('M') = 128;
        v_('N') = 3;
        v_('1') = 1;
        v_('X11') = 1.0;
        v_('X12') = 1.0;
        v_('X13') = 1.0;
        v_('X14') = 1.0;
        v_('X15') = 1.0;
        v_('X16') = 1.0;
        v_('X17') = 1.0;
        v_('X18') = 1.0;
        v_('X19') = 1.0;
        v_('X110') = 1.0;
        v_('X111') = 1.0;
        v_('X112') = 1.0;
        v_('X113') = 1.0;
        v_('X114') = 1.0;
        v_('X115') = 1.0;
        v_('X116') = 1.0;
        v_('X117') = 2.0;
        v_('X118') = 2.0;
        v_('X119') = 2.0;
        v_('X120') = 2.0;
        v_('X121') = 2.0;
        v_('X122') = 2.0;
        v_('X123') = 2.0;
        v_('X124') = 2.0;
        v_('X125') = 2.0;
        v_('X126') = 2.0;
        v_('X127') = 2.0;
        v_('X128') = 2.0;
        v_('X129') = 2.0;
        v_('X130') = 2.0;
        v_('X131') = 2.0;
        v_('X132') = 2.0;
        v_('X133') = 4.0;
        v_('X134') = 4.0;
        v_('X135') = 4.0;
        v_('X136') = 4.0;
        v_('X137') = 4.0;
        v_('X138') = 4.0;
        v_('X139') = 4.0;
        v_('X140') = 4.0;
        v_('X141') = 4.0;
        v_('X142') = 4.0;
        v_('X143') = 4.0;
        v_('X144') = 4.0;
        v_('X145') = 4.0;
        v_('X146') = 4.0;
        v_('X147') = 4.0;
        v_('X148') = 4.0;
        v_('X149') = 8.0;
        v_('X150') = 8.0;
        v_('X151') = 8.0;
        v_('X152') = 8.0;
        v_('X153') = 8.0;
        v_('X154') = 8.0;
        v_('X155') = 8.0;
        v_('X156') = 8.0;
        v_('X157') = 8.0;
        v_('X158') = 8.0;
        v_('X159') = 8.0;
        v_('X160') = 8.0;
        v_('X161') = 8.0;
        v_('X162') = 8.0;
        v_('X163') = 8.0;
        v_('X164') = 8.0;
        v_('X165') = 16.0;
        v_('X166') = 16.0;
        v_('X167') = 16.0;
        v_('X168') = 16.0;
        v_('X169') = 16.0;
        v_('X170') = 16.0;
        v_('X171') = 16.0;
        v_('X172') = 16.0;
        v_('X173') = 16.0;
        v_('X174') = 16.0;
        v_('X175') = 16.0;
        v_('X176') = 16.0;
        v_('X177') = 16.0;
        v_('X178') = 16.0;
        v_('X179') = 16.0;
        v_('X180') = 16.0;
        v_('X181') = 32.0;
        v_('X182') = 32.0;
        v_('X183') = 32.0;
        v_('X184') = 32.0;
        v_('X185') = 32.0;
        v_('X186') = 32.0;
        v_('X187') = 32.0;
        v_('X188') = 32.0;
        v_('X189') = 32.0;
        v_('X190') = 32.0;
        v_('X191') = 32.0;
        v_('X192') = 32.0;
        v_('X193') = 32.0;
        v_('X194') = 32.0;
        v_('X195') = 32.0;
        v_('X196') = 32.0;
        v_('X197') = 48.0;
        v_('X198') = 48.0;
        v_('X199') = 48.0;
        v_('X1100') = 48.0;
        v_('X1101') = 48.0;
        v_('X1102') = 48.0;
        v_('X1103') = 48.0;
        v_('X1104') = 48.0;
        v_('X1105') = 48.0;
        v_('X1106') = 48.0;
        v_('X1107') = 48.0;
        v_('X1108') = 48.0;
        v_('X1109') = 48.0;
        v_('X1110') = 48.0;
        v_('X1111') = 48.0;
        v_('X1112') = 48.0;
        v_('X1113') = 64.0;
        v_('X1114') = 64.0;
        v_('X1115') = 64.0;
        v_('X1116') = 64.0;
        v_('X1117') = 64.0;
        v_('X1118') = 64.0;
        v_('X1119') = 64.0;
        v_('X1120') = 64.0;
        v_('X1121') = 64.0;
        v_('X1122') = 64.0;
        v_('X1123') = 64.0;
        v_('X1124') = 64.0;
        v_('X1125') = 64.0;
        v_('X1126') = 64.0;
        v_('X1127') = 64.0;
        v_('X1128') = 64.0;
        v_('X21') = 180.0;
        v_('X22') = 180.0;
        v_('X23') = 180.0;
        v_('X24') = 180.0;
        v_('X25') = 225.0;
        v_('X26') = 225.0;
        v_('X27') = 225.0;
        v_('X28') = 225.0;
        v_('X29') = 250.0;
        v_('X210') = 250.0;
        v_('X211') = 250.0;
        v_('X212') = 250.0;
        v_('X213') = 275.0;
        v_('X214') = 275.0;
        v_('X215') = 275.0;
        v_('X216') = 275.0;
        v_('X217') = 180.0;
        v_('X218') = 180.0;
        v_('X219') = 180.0;
        v_('X220') = 180.0;
        v_('X221') = 225.0;
        v_('X222') = 225.0;
        v_('X223') = 225.0;
        v_('X224') = 225.0;
        v_('X225') = 250.0;
        v_('X226') = 250.0;
        v_('X227') = 250.0;
        v_('X228') = 250.0;
        v_('X229') = 275.0;
        v_('X230') = 275.0;
        v_('X231') = 275.0;
        v_('X232') = 275.0;
        v_('X233') = 180.0;
        v_('X234') = 180.0;
        v_('X235') = 180.0;
        v_('X236') = 180.0;
        v_('X237') = 225.0;
        v_('X238') = 225.0;
        v_('X239') = 225.0;
        v_('X240') = 225.0;
        v_('X241') = 250.0;
        v_('X242') = 250.0;
        v_('X243') = 250.0;
        v_('X244') = 250.0;
        v_('X245') = 275.0;
        v_('X246') = 275.0;
        v_('X247') = 275.0;
        v_('X248') = 275.0;
        v_('X249') = 180.0;
        v_('X250') = 180.0;
        v_('X251') = 180.0;
        v_('X252') = 180.0;
        v_('X253') = 225.0;
        v_('X254') = 225.0;
        v_('X255') = 225.0;
        v_('X256') = 225.0;
        v_('X257') = 250.0;
        v_('X258') = 250.0;
        v_('X259') = 250.0;
        v_('X260') = 250.0;
        v_('X261') = 275.0;
        v_('X262') = 275.0;
        v_('X263') = 275.0;
        v_('X264') = 275.0;
        v_('X265') = 180.0;
        v_('X266') = 180.0;
        v_('X267') = 180.0;
        v_('X268') = 180.0;
        v_('X269') = 225.0;
        v_('X270') = 225.0;
        v_('X271') = 225.0;
        v_('X272') = 225.0;
        v_('X273') = 250.0;
        v_('X274') = 250.0;
        v_('X275') = 250.0;
        v_('X276') = 250.0;
        v_('X277') = 275.0;
        v_('X278') = 275.0;
        v_('X279') = 275.0;
        v_('X280') = 275.0;
        v_('X281') = 180.0;
        v_('X282') = 180.0;
        v_('X283') = 180.0;
        v_('X284') = 180.0;
        v_('X285') = 225.0;
        v_('X286') = 225.0;
        v_('X287') = 225.0;
        v_('X288') = 225.0;
        v_('X289') = 250.0;
        v_('X290') = 250.0;
        v_('X291') = 250.0;
        v_('X292') = 250.0;
        v_('X293') = 275.0;
        v_('X294') = 275.0;
        v_('X295') = 275.0;
        v_('X296') = 275.0;
        v_('X297') = 180.0;
        v_('X298') = 180.0;
        v_('X299') = 180.0;
        v_('X2100') = 180.0;
        v_('X2101') = 225.0;
        v_('X2102') = 225.0;
        v_('X2103') = 225.0;
        v_('X2104') = 225.0;
        v_('X2105') = 250.0;
        v_('X2106') = 250.0;
        v_('X2107') = 250.0;
        v_('X2108') = 250.0;
        v_('X2109') = 275.0;
        v_('X2110') = 275.0;
        v_('X2111') = 275.0;
        v_('X2112') = 275.0;
        v_('X2113') = 180.0;
        v_('X2114') = 180.0;
        v_('X2115') = 180.0;
        v_('X2116') = 180.0;
        v_('X2117') = 225.0;
        v_('X2118') = 225.0;
        v_('X2119') = 225.0;
        v_('X2120') = 225.0;
        v_('X2121') = 250.0;
        v_('X2122') = 250.0;
        v_('X2123') = 250.0;
        v_('X2124') = 250.0;
        v_('X2125') = 275.0;
        v_('X2126') = 275.0;
        v_('X2127') = 275.0;
        v_('X2128') = 275.0;
        v_('Y1') = 15.00;
        v_('Y2') = 17.00;
        v_('Y3') = 15.50;
        v_('Y4') = 16.50;
        v_('Y5') = 15.50;
        v_('Y6') = 15.00;
        v_('Y7') = 16.00;
        v_('Y8') = 14.50;
        v_('Y9') = 15.00;
        v_('Y10') = 14.50;
        v_('Y11') = 12.50;
        v_('Y12') = 11.00;
        v_('Y13') = 14.00;
        v_('Y14') = 13.00;
        v_('Y15') = 14.00;
        v_('Y16') = 11.50;
        v_('Y17') = 14.00;
        v_('Y18') = 16.00;
        v_('Y19') = 13.00;
        v_('Y20') = 13.50;
        v_('Y21') = 13.00;
        v_('Y22') = 13.50;
        v_('Y23') = 12.50;
        v_('Y24') = 12.50;
        v_('Y25') = 12.50;
        v_('Y26') = 12.00;
        v_('Y27') = 11.50;
        v_('Y28') = 12.00;
        v_('Y29') = 13.00;
        v_('Y30') = 11.50;
        v_('Y31') = 13.00;
        v_('Y32') = 12.50;
        v_('Y33') = 13.50;
        v_('Y34') = 17.50;
        v_('Y35') = 17.50;
        v_('Y36') = 13.50;
        v_('Y37') = 12.50;
        v_('Y38') = 12.50;
        v_('Y39') = 15.00;
        v_('Y40') = 13.00;
        v_('Y41') = 12.00;
        v_('Y42') = 13.00;
        v_('Y43') = 12.00;
        v_('Y44') = 13.50;
        v_('Y45') = 10.00;
        v_('Y46') = 11.50;
        v_('Y47') = 11.00;
        v_('Y48') = 9.50;
        v_('Y49') = 15.00;
        v_('Y50') = 15.00;
        v_('Y51') = 15.50;
        v_('Y52') = 16.00;
        v_('Y53') = 13.00;
        v_('Y54') = 10.50;
        v_('Y55') = 13.50;
        v_('Y56') = 14.00;
        v_('Y57') = 12.50;
        v_('Y58') = 12.00;
        v_('Y59') = 11.50;
        v_('Y60') = 11.50;
        v_('Y61') = 6.50;
        v_('Y62') = 5.50;
        v_('Y63') = 6.00;
        v_('Y64') = 6.00;
        v_('Y65') = 18.50;
        v_('Y66') = 17.00;
        v_('Y67') = 15.30;
        v_('Y68') = 16.00;
        v_('Y69') = 13.00;
        v_('Y70') = 14.00;
        v_('Y71') = 12.50;
        v_('Y72') = 11.00;
        v_('Y73') = 12.00;
        v_('Y74') = 12.00;
        v_('Y75') = 11.50;
        v_('Y76') = 12.00;
        v_('Y77') = 6.00;
        v_('Y78') = 6.00;
        v_('Y79') = 5.00;
        v_('Y80') = 5.50;
        v_('Y81') = 12.50;
        v_('Y82') = 13.00;
        v_('Y83') = 16.00;
        v_('Y84') = 12.00;
        v_('Y85') = 11.00;
        v_('Y86') = 9.50;
        v_('Y87') = 11.00;
        v_('Y88') = 11.00;
        v_('Y89') = 11.00;
        v_('Y90') = 10.00;
        v_('Y91') = 10.50;
        v_('Y92') = 10.50;
        v_('Y93') = 2.70;
        v_('Y94') = 2.70;
        v_('Y95') = 2.50;
        v_('Y96') = 2.40;
        v_('Y97') = 13.00;
        v_('Y98') = 13.50;
        v_('Y99') = 16.50;
        v_('Y100') = 13.60;
        v_('Y101') = 11.50;
        v_('Y102') = 10.50;
        v_('Y103') = 13.50;
        v_('Y104') = 12.00;
        v_('Y105') = 7.00;
        v_('Y106') = 6.90;
        v_('Y107') = 8.80;
        v_('Y108') = 7.90;
        v_('Y109') = 1.20;
        v_('Y110') = 1.50;
        v_('Y111') = 1.00;
        v_('Y112') = 1.50;
        v_('Y113') = 13.00;
        v_('Y114') = 12.50;
        v_('Y115') = 16.50;
        v_('Y116') = 16.00;
        v_('Y117') = 11.00;
        v_('Y118') = 11.50;
        v_('Y119') = 10.50;
        v_('Y120') = 10.00;
        v_('Y121') = 7.27;
        v_('Y122') = 7.50;
        v_('Y123') = 6.70;
        v_('Y124') = 7.60;
        v_('Y125') = 1.50;
        v_('Y126') = 1.00;
        v_('Y127') = 1.20;
        v_('Y128') = 1.20;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['F',int2str(I)];
            iv = ix_('B1');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
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
        for I=v_('1'):v_('M')
            v_('LOGY') = log(v_(['Y',int2str(I)]));
            pbm.gconst(ig_(['F',int2str(I)])) = v_('LOGY');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'B1'))
            pb.x0(ix_('B1'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('B1')),1) = 2.0;
        end
        if(isKey(ix_,'B2'))
            pb.x0(ix_('B2'),1) = 0.0001;
        else
            pb.y0(find(pbm.congrps==ig_('B2')),1) = 0.0001;
        end
        if(isKey(ix_,'B3'))
            pb.x0(ix_('B3'),1) = -0.01;
        else
            pb.y0(find(pbm.congrps==ig_('B3')),1) = -0.01;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE6',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'X1';
        elftp{it}{2} = 'X2';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE6';
            ielftype(ie) = iet_('eE6');
            vname = 'B2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X1',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X1',int2str(I)]);
            [~,posep] = ismember('X2',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X2',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-MN-3-128';
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

    case 'eE6'

        EV_  = varargin{1};
        iel_ = varargin{2};
        E = exp(-EV_(2)*pbm.elpar{iel_}(2));
        X1E = pbm.elpar{iel_}(1)*E;
        V1X1E = EV_(1)*pbm.elpar{iel_}(1)*E;
        varargout{1} = V1X1E;
        if(nargout>1)
            g_(1,1) = X1E;
            g_(2,1) = -V1X1E*pbm.elpar{iel_}(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -X1E*pbm.elpar{iel_}(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = V1X1E*pbm.elpar{iel_}(2)^2;
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

