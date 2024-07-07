function varargout = MAXLIKA(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MAXLIKA
%    *********
% 
%    A variant of Hock and Schittkowski problem 105, where the
%    (inactive) inequality constraint is dropped.
% 
%    Source:
%    Ph. Toint and A. Griewank.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'OBR2-AY-8-0'
% 
%    Constants
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MAXLIKA';

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
        v_('1') = 1;
        v_('235') = 235;
        v_('Y1') = 95.0;
        v_('Y2') = 105.0;
        v_('Y3') = 110.0;
        v_('Y4') = 110.0;
        v_('Y5') = 110.0;
        v_('Y6') = 110.0;
        v_('Y7') = 115.0;
        v_('Y8') = 115.0;
        v_('Y9') = 115.0;
        v_('Y10') = 115.0;
        v_('Y11') = 120.0;
        v_('Y12') = 120.0;
        v_('Y13') = 120.0;
        v_('Y14') = 120.0;
        v_('Y15') = 120.0;
        v_('Y16') = 120.0;
        v_('Y17') = 120.0;
        v_('Y18') = 120.0;
        v_('Y19') = 120.0;
        v_('Y20') = 120.0;
        v_('Y21') = 120.0;
        v_('Y22') = 120.0;
        v_('Y23') = 120.0;
        v_('Y24') = 120.0;
        v_('Y25') = 120.0;
        v_('Y26') = 125.0;
        v_('Y27') = 125.0;
        v_('Y28') = 125.0;
        v_('Y29') = 125.0;
        v_('Y30') = 125.0;
        v_('Y31') = 125.0;
        v_('Y32') = 125.0;
        v_('Y33') = 125.0;
        v_('Y34') = 125.0;
        v_('Y35') = 125.0;
        v_('Y36') = 125.0;
        v_('Y37') = 125.0;
        v_('Y38') = 125.0;
        v_('Y39') = 125.0;
        v_('Y40') = 125.0;
        v_('Y41') = 130.0;
        v_('Y42') = 130.0;
        v_('Y43') = 130.0;
        v_('Y44') = 130.0;
        v_('Y45') = 130.0;
        v_('Y46') = 130.0;
        v_('Y47') = 130.0;
        v_('Y48') = 130.0;
        v_('Y49') = 130.0;
        v_('Y50') = 130.0;
        v_('Y51') = 130.0;
        v_('Y52') = 130.0;
        v_('Y53') = 130.0;
        v_('Y54') = 130.0;
        v_('Y55') = 130.0;
        v_('Y56') = 135.0;
        v_('Y57') = 135.0;
        v_('Y58') = 135.0;
        v_('Y59') = 135.0;
        v_('Y60') = 135.0;
        v_('Y61') = 135.0;
        v_('Y62') = 135.0;
        v_('Y63') = 135.0;
        v_('Y64') = 135.0;
        v_('Y65') = 135.0;
        v_('Y66') = 135.0;
        v_('Y67') = 135.0;
        v_('Y68') = 135.0;
        v_('Y69') = 140.0;
        v_('Y70') = 140.0;
        v_('Y71') = 140.0;
        v_('Y72') = 140.0;
        v_('Y73') = 140.0;
        v_('Y74') = 140.0;
        v_('Y75') = 140.0;
        v_('Y76') = 140.0;
        v_('Y77') = 140.0;
        v_('Y78') = 140.0;
        v_('Y79') = 140.0;
        v_('Y80') = 140.0;
        v_('Y81') = 140.0;
        v_('Y82') = 140.0;
        v_('Y83') = 140.0;
        v_('Y84') = 140.0;
        v_('Y85') = 140.0;
        v_('Y86') = 140.0;
        v_('Y87') = 140.0;
        v_('Y88') = 140.0;
        v_('Y89') = 140.0;
        v_('Y90') = 145.0;
        v_('Y91') = 145.0;
        v_('Y92') = 145.0;
        v_('Y93') = 145.0;
        v_('Y94') = 145.0;
        v_('Y95') = 145.0;
        v_('Y96') = 145.0;
        v_('Y97') = 145.0;
        v_('Y98') = 145.0;
        v_('Y99') = 145.0;
        v_('Y100') = 145.0;
        v_('Y101') = 145.0;
        v_('Y102') = 150.0;
        v_('Y103') = 150.0;
        v_('Y104') = 150.0;
        v_('Y105') = 150.0;
        v_('Y106') = 150.0;
        v_('Y107') = 150.0;
        v_('Y108') = 150.0;
        v_('Y109') = 150.0;
        v_('Y110') = 150.0;
        v_('Y111') = 150.0;
        v_('Y112') = 150.0;
        v_('Y113') = 150.0;
        v_('Y114') = 150.0;
        v_('Y115') = 150.0;
        v_('Y116') = 150.0;
        v_('Y117') = 150.0;
        v_('Y118') = 150.0;
        v_('Y119') = 155.0;
        v_('Y120') = 155.0;
        v_('Y121') = 155.0;
        v_('Y122') = 155.0;
        v_('Y123') = 160.0;
        v_('Y124') = 160.0;
        v_('Y125') = 160.0;
        v_('Y126') = 160.0;
        v_('Y127') = 160.0;
        v_('Y128') = 160.0;
        v_('Y129') = 160.0;
        v_('Y130') = 160.0;
        v_('Y131') = 160.0;
        v_('Y132') = 160.0;
        v_('Y133') = 160.0;
        v_('Y134') = 160.0;
        v_('Y135') = 160.0;
        v_('Y136') = 160.0;
        v_('Y137') = 160.0;
        v_('Y138') = 160.0;
        v_('Y139') = 160.0;
        v_('Y140') = 160.0;
        v_('Y141') = 160.0;
        v_('Y142') = 160.0;
        v_('Y143') = 165.0;
        v_('Y144') = 165.0;
        v_('Y145') = 165.0;
        v_('Y146') = 165.0;
        v_('Y147') = 165.0;
        v_('Y148') = 165.0;
        v_('Y149') = 165.0;
        v_('Y150') = 165.0;
        v_('Y151') = 170.0;
        v_('Y152') = 170.0;
        v_('Y153') = 170.0;
        v_('Y154') = 170.0;
        v_('Y155') = 170.0;
        v_('Y156') = 170.0;
        v_('Y157') = 170.0;
        v_('Y158') = 170.0;
        v_('Y159') = 170.0;
        v_('Y160') = 170.0;
        v_('Y161') = 170.0;
        v_('Y162') = 170.0;
        v_('Y163') = 170.0;
        v_('Y164') = 170.0;
        v_('Y165') = 170.0;
        v_('Y166') = 170.0;
        v_('Y167') = 170.0;
        v_('Y168') = 175.0;
        v_('Y169') = 175.0;
        v_('Y170') = 175.0;
        v_('Y171') = 175.0;
        v_('Y172') = 175.0;
        v_('Y173') = 175.0;
        v_('Y174') = 175.0;
        v_('Y175') = 175.0;
        v_('Y176') = 180.0;
        v_('Y177') = 180.0;
        v_('Y178') = 180.0;
        v_('Y179') = 180.0;
        v_('Y180') = 180.0;
        v_('Y181') = 180.0;
        v_('Y182') = 185.0;
        v_('Y183') = 185.0;
        v_('Y184') = 185.0;
        v_('Y185') = 185.0;
        v_('Y186') = 185.0;
        v_('Y187') = 185.0;
        v_('Y188') = 190.0;
        v_('Y189') = 190.0;
        v_('Y190') = 190.0;
        v_('Y191') = 190.0;
        v_('Y192') = 190.0;
        v_('Y193') = 190.0;
        v_('Y194') = 190.0;
        v_('Y195') = 195.0;
        v_('Y196') = 195.0;
        v_('Y197') = 195.0;
        v_('Y198') = 195.0;
        v_('Y199') = 200.0;
        v_('Y200') = 200.0;
        v_('Y201') = 200.0;
        v_('Y202') = 205.0;
        v_('Y203') = 205.0;
        v_('Y204') = 205.0;
        v_('Y205') = 210.0;
        v_('Y206') = 210.0;
        v_('Y207') = 210.0;
        v_('Y208') = 210.0;
        v_('Y209') = 210.0;
        v_('Y210') = 210.0;
        v_('Y211') = 210.0;
        v_('Y212') = 210.0;
        v_('Y213') = 215.0;
        v_('Y214') = 220.0;
        v_('Y215') = 220.0;
        v_('Y216') = 220.0;
        v_('Y217') = 220.0;
        v_('Y218') = 220.0;
        v_('Y219') = 220.0;
        v_('Y220') = 230.0;
        v_('Y221') = 230.0;
        v_('Y222') = 230.0;
        v_('Y223') = 230.0;
        v_('Y224') = 230.0;
        v_('Y225') = 235.0;
        v_('Y226') = 240.0;
        v_('Y227') = 240.0;
        v_('Y228') = 240.0;
        v_('Y229') = 240.0;
        v_('Y230') = 240.0;
        v_('Y231') = 240.0;
        v_('Y232') = 240.0;
        v_('Y233') = 245.0;
        v_('Y234') = 250.0;
        v_('Y235') = 250.0;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','X1',ix_);
        pb.xnames{iv} = 'X1';
        [iv,ix_] = s2mpjlib('ii','X2',ix_);
        pb.xnames{iv} = 'X2';
        [iv,ix_] = s2mpjlib('ii','X3',ix_);
        pb.xnames{iv} = 'X3';
        [iv,ix_] = s2mpjlib('ii','X4',ix_);
        pb.xnames{iv} = 'X4';
        [iv,ix_] = s2mpjlib('ii','X5',ix_);
        pb.xnames{iv} = 'X5';
        [iv,ix_] = s2mpjlib('ii','X6',ix_);
        pb.xnames{iv} = 'X6';
        [iv,ix_] = s2mpjlib('ii','X7',ix_);
        pb.xnames{iv} = 'X7';
        [iv,ix_] = s2mpjlib('ii','X8',ix_);
        pb.xnames{iv} = 'X8';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('235')
            [ig,ig_] = s2mpjlib('ii',['L',int2str(I)],ig_);
            gtype{ig} = '<>';
            pbm.gscale(ig,1) = -1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('X1'),1) = 0.001;
        pb.xupper(ix_('X1')) = 0.499;
        pb.xlower(ix_('X2'),1) = 0.001;
        pb.xupper(ix_('X2')) = 0.499;
        pb.xlower(ix_('X3'),1) = 100.0;
        pb.xupper(ix_('X3')) = 180.0;
        pb.xlower(ix_('X4'),1) = 130.0;
        pb.xupper(ix_('X4')) = 210.0;
        pb.xlower(ix_('X5'),1) = 170.0;
        pb.xupper(ix_('X5')) = 240.0;
        pb.xlower(ix_('X6'),1) = 5.0;
        pb.xupper(ix_('X6')) = 25.0;
        pb.xlower(ix_('X7'),1) = 5.0;
        pb.xupper(ix_('X7')) = 25.0;
        pb.xlower(ix_('X8'),1) = 5.0;
        pb.xupper(ix_('X8')) = 25.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('X1'),1) = 0.1;
        pb.x0(ix_('X2'),1) = 0.2;
        pb.x0(ix_('X3'),1) = 100.0;
        pb.x0(ix_('X4'),1) = 125.0;
        pb.x0(ix_('X5'),1) = 175.0;
        pb.x0(ix_('X6'),1) = 11.2;
        pb.x0(ix_('X7'),1) = 13.2;
        pb.x0(ix_('X8'),1) = 15.8;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eAB',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'V';
        elftv{it}{3} = 'W';
        elftp{it}{1} = 'Y';
        [it,iet_] = s2mpjlib( 'ii', 'eC',iet_);
        elftv{it}{1} = 'U';
        elftv{it}{2} = 'X';
        elftv{it}{3} = 'S';
        elftv{it}{4} = 'T';
        elftp{it}{1} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('235')
            ename = ['A',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eAB';
            ielftype(ie) = iet_('eAB');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('Y',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
            ename = ['B',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eAB';
            ielftype(ie) = iet_('eAB');
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X7';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('Y',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
            ename = ['C',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eC';
            ielftype(ie) = iet_('eC');
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('U',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X8';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('S',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('T',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('Y',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gLN',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for I=v_('1'):v_('235')
            ig = ig_(['L',int2str(I)]);
            pbm.grftype{ig} = 'gLN';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['A',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['B',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['C',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'OBR2-AY-8-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eAB'

        EV_  = varargin{1};
        iel_ = varargin{2};
        YMW = pbm.elpar{iel_}(1)-EV_(3);
        YMWSQ = YMW*YMW;
        VSQ = EV_(2)*EV_(2);
        VCB = VSQ*EV_(2);
        A = -YMWSQ/(2.0*VSQ);
        DADV = YMWSQ/VCB;
        DADW = YMW/VSQ;
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ);
        D2ADVW = -2.0*YMW/VCB;
        D2ADW2 = -1.0/VSQ;
        E = exp(A);
        DEDV = E*DADV;
        DEDW = E*DADW;
        B = EV_(1)*E;
        DBDV = B*DADV;
        DBDW = B*DADW;
        D2BDV2 = DBDV*DADV+B*D2ADV2;
        D2BDVW = DBDW*DADV+B*D2ADVW;
        D2BDW2 = DBDW*DADW+B*D2ADW2;
        varargout{1} = B/EV_(2);
        if(nargout>1)
            g_(1,1) = E/EV_(2);
            g_(2,1) = (DBDV-B/EV_(2))/EV_(2);
            g_(3,1) = DBDW/EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = (DEDV-E/EV_(2))/EV_(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = DEDW/EV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = (D2BDV2-DBDV/EV_(2)+B/VSQ)/EV_(2)-(DBDV-B/EV_(2))/VSQ;
                H_(2,3) = (D2BDVW-DBDW/EV_(2))/EV_(2);
                H_(3,2) = H_(2,3);
                H_(3,3) = D2BDW2/EV_(2);
                varargout{3} = H_;
            end
        end

    case 'eC'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(3,4);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)+1;
        U_(2,3) = U_(2,3)+1;
        U_(3,4) = U_(3,4)+1;
        IV_(1) = U_(1,:)*EV_;
        IV_(2) = U_(2,:)*EV_;
        IV_(3) = U_(3,:)*EV_;
        YMW = pbm.elpar{iel_}(1)-IV_(3);
        YMWSQ = YMW*YMW;
        VSQ = IV_(2)*IV_(2);
        VCB = VSQ*IV_(2);
        A = -YMWSQ/(2.0*VSQ);
        DADV = YMWSQ/VCB;
        DADW = YMW/VSQ;
        D2ADV2 = -3.0*YMWSQ/(VSQ*VSQ);
        D2ADVW = -2.0*YMW/VCB;
        D2ADW2 = -1.0/VSQ;
        E = exp(A);
        DEDV = E*DADV;
        DEDW = E*DADW;
        B = (1.0-IV_(1))*E;
        DBDZ = -E;
        DBDV = B*DADV;
        DBDW = B*DADW;
        D2BDVZ = -DEDV;
        D2BDWZ = -DEDW;
        D2BDV2 = DBDV*DADV+B*D2ADV2;
        D2BDVW = DBDW*DADV+B*D2ADVW;
        D2BDW2 = DBDW*DADW+B*D2ADW2;
        varargout{1} = B/IV_(2);
        if(nargout>1)
            g_(1,1) = DBDZ/IV_(2);
            g_(2,1) = (DBDV-B/IV_(2))/IV_(2);
            g_(3,1) = DBDW/IV_(2);
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = (DEDV-E/IV_(2))/IV_(2);
                H_(2,1) = H_(1,2);
                H_(1,3) = D2BDWZ/IV_(2);
                H_(3,1) = H_(1,3);
                H_(2,2) = (D2BDV2-DBDV/IV_(2)+B/VSQ)/IV_(2)-(DBDV-B/IV_(2))/VSQ;
                H_(2,3) = (D2BDVW-DBDW/IV_(2))/IV_(2);
                H_(3,2) = H_(2,3);
                H_(3,3) = D2BDW2/IV_(2);
                varargout{3} = U_.'*H_*U_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gLN'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = log(GVAR_*0.39894228);
        if(nargout>1)
            g_ = 1.0/GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = -1.0/GVAR_^2;
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

