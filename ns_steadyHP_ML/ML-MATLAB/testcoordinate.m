function testcoordinate()
v_coord =[

   -0.0456    0.5008
    0.0017    0.4259
    0.1136    0.5253
    0.1750    0.4122
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0
         0         0];


nEv =9;


nR =3;
get_cordinates(v_coord,nEv,nR)

end

function v_coord = get_cordinates(v_coord,nEv,nR)
%nEv
tmpindex = [1,3;2,4;1,2;3,4];
if nEv>2
    for edge = 1:4
        evstart = 5+(nEv-2)*(edge-1);
        %evstart:evstart+nEv-3
        cc = 1/(nEv-1);
        count = 1;
        index1 = tmpindex(edge,1);
        index2 = tmpindex(edge,2);
        for vv = evstart:evstart+nEv-3
            c2 = count*cc;
            c1 = 1-c2;
            v_coord(vv,:) = c1*v_coord(index1,:)+c2*v_coord(index2,:);
            %eV(edge,3:nEv) = vID(evstart:evstart+nEv-3);
            count = count+1;
        end
    end
    
    %nFace = (nEv-2)^2;
    if nR == 1
        v_coord(9,:) = 0.25*(v_coord(5,:)+v_coord(6,:)+v_coord(7,:)+v_coord(8,:));
    end
    
    if nR ==2
        v_coord(21,:) = 0.25*(v_coord(6,:)+v_coord(9,:)+v_coord(12,:)+v_coord(15,:));
        
        v_coord(18,:) = 0.5*(v_coord(12,:)+v_coord(21,:));
        v_coord(20,:) = 0.5*(v_coord(6,:)+v_coord(21,:));
        v_coord(22,:) = 0.5*(v_coord(21,:)+v_coord(9,:));
        v_coord(24,:) = 0.5*(v_coord(21,:)+v_coord(15,:));
        
        v_coord(17,:) = 0.25*(v_coord(5,:)+v_coord(18,:)+v_coord(11,:)+v_coord(20,:));
        v_coord(19,:) = 0.25*(v_coord(18,:)+v_coord(8,:)+v_coord(13,:)+v_coord(22,:));
        v_coord(23,:) = 0.25*(v_coord(20,:)+v_coord(14,:)+v_coord(7,:)+v_coord(24,:));
        v_coord(25,:) = 0.25*(v_coord(22,:)+v_coord(16,:)+v_coord(24,:)+v_coord(10,:));

    end
    
    if nR ==3
        v_coord(57,:) = 0.25*(v_coord(8,:)+v_coord(15,:)+v_coord(22,:)+v_coord(29,:));
        
        v_coord(55,:) = 0.5*(v_coord(8,:)+v_coord(57,:));
        v_coord(59,:) = 0.5*(v_coord(15,:)+v_coord(57,:));
        v_coord(43,:) = 0.5*(v_coord(22,:)+v_coord(57,:));
        v_coord(71,:) = 0.5*(v_coord(29,:)+v_coord(57,:));
        
        v_coord(41,:) = 0.25*(v_coord(6,:)+v_coord(43,:)+v_coord(20,:)+v_coord(55,:));
        v_coord(45,:) = 0.25*(v_coord(43,:)+v_coord(13,:)+v_coord(24,:)+v_coord(59,:));
        v_coord(69,:) = 0.25*(v_coord(10,:)+v_coord(71,:)+v_coord(55,:)+v_coord(27,:));
        v_coord(73,:) = 0.25*(v_coord(71,:)+v_coord(17,:)+v_coord(59,:)+v_coord(31,:));
        
         v_coord(34,:) = 0.5*(v_coord(20,:)+v_coord(41,:));
         v_coord(36,:) = 0.5*(v_coord(22,:)+v_coord(43,:));
         v_coord(38,:) = 0.5*(v_coord(24,:)+v_coord(45,:));  
         
         v_coord(40,:) = 0.5*(v_coord(6,:)+v_coord(41,:));
         v_coord(42,:) = 0.5*(v_coord(41,:)+v_coord(43,:));
         v_coord(44,:) = 0.5*(v_coord(43,:)+v_coord(45,:));
         v_coord(46,:) = 0.5*(v_coord(45,:)+v_coord(13,:));
         
         v_coord(48,:) = 0.5*(v_coord(41,:)+v_coord(55,:));
         v_coord(50,:) = 0.5*(v_coord(43,:)+v_coord(57,:));
         v_coord(52,:) = 0.5*(v_coord(45,:)+v_coord(59,:));  
        
         v_coord(54,:) = 0.5*(v_coord(8,:)+v_coord(55,:));
         v_coord(56,:) = 0.5*(v_coord(55,:)+v_coord(57,:));
         v_coord(58,:) = 0.5*(v_coord(57,:)+v_coord(59,:));
         v_coord(60,:) = 0.5*(v_coord(59,:)+v_coord(15,:));
         
         v_coord(62,:) = 0.5*(v_coord(55,:)+v_coord(69,:));
         v_coord(64,:) = 0.5*(v_coord(57,:)+v_coord(71,:));
         v_coord(66,:) = 0.5*(v_coord(59,:)+v_coord(73,:));  
        
         v_coord(68,:) = 0.5*(v_coord(10,:)+v_coord(69,:));       
         v_coord(70,:) = 0.5*(v_coord(69,:)+v_coord(71,:));
         v_coord(72,:) = 0.5*(v_coord(71,:)+v_coord(73,:));
         v_coord(74,:) = 0.5*(v_coord(73,:)+v_coord(17,:));
         
         v_coord(76,:) = 0.5*(v_coord(69,:)+v_coord(27,:));
         v_coord(78,:) = 0.5*(v_coord(71,:)+v_coord(29,:));
         v_coord(80,:) = 0.5*(v_coord(73,:)+v_coord(31,:));  
        
        v_coord(33,:) = 0.25*(v_coord(19,:)+v_coord(40,:)+v_coord(5,:)+v_coord(34,:));
        v_coord(35,:) = 0.25*(v_coord(21,:)+v_coord(42,:)+v_coord(34,:)+v_coord(36,:));
        v_coord(37,:) = 0.25*(v_coord(23,:)+v_coord(44,:)+v_coord(36,:)+v_coord(38,:));
        v_coord(39,:) = 0.25*(v_coord(25,:)+v_coord(46,:)+v_coord(38,:)+v_coord(12,:));

        v_coord(47,:) = 0.25*(v_coord(7,:)+v_coord(48,:)+v_coord(40,:)+v_coord(54,:));
        v_coord(49,:) = 0.25*(v_coord(48,:)+v_coord(50,:)+v_coord(42,:)+v_coord(56,:));
        v_coord(51,:) = 0.25*(v_coord(50,:)+v_coord(52,:)+v_coord(44,:)+v_coord(58,:));
        v_coord(53,:) = 0.25*(v_coord(52,:)+v_coord(14,:)+v_coord(46,:)+v_coord(60,:));
        
        v_coord(61,:) = 0.25*(v_coord(9,:)+v_coord(62,:)+v_coord(54,:)+v_coord(68,:));
        v_coord(63,:) = 0.25*(v_coord(70,:)+v_coord(56,:)+v_coord(62,:)+v_coord(64,:));
        v_coord(65,:) = 0.25*(v_coord(58,:)+v_coord(72,:)+v_coord(64,:)+v_coord(66,:));
        v_coord(67,:) = 0.25*(v_coord(60,:)+v_coord(74,:)+v_coord(66,:)+v_coord(16,:));
        
        v_coord(75,:) = 0.25*(v_coord(68,:)+v_coord(26,:)+v_coord(11,:)+v_coord(76,:));
        v_coord(77,:) = 0.25*(v_coord(28,:)+v_coord(70,:)+v_coord(76,:)+v_coord(78,:));
        v_coord(79,:) = 0.25*(v_coord(30,:)+v_coord(72,:)+v_coord(78,:)+v_coord(80,:));
        v_coord(81,:) = 0.25*(v_coord(32,:)+v_coord(74,:)+v_coord(80,:)+v_coord(18,:));
    end
    
    if nR ==4
        v_coord(177,:) = 0.25*(v_coord(12,:)+v_coord(27,:)+v_coord(42,:)+v_coord(57,:));
        
        v_coord(173,:) = 0.5*(v_coord(12,:)+v_coord(177,:));
        v_coord(181,:) = 0.5*(v_coord(27,:)+v_coord(177,:));
        v_coord(117,:) = 0.5*(v_coord(42,:)+v_coord(177,:));
        v_coord(237,:) = 0.5*(v_coord(57,:)+v_coord(177,:));
        
        v_coord(113,:) = 0.25*(v_coord(8,:)+v_coord(117,:)+v_coord(38,:)+v_coord(173,:));
        v_coord(121,:) = 0.25*(v_coord(117,:)+v_coord(23,:)+v_coord(46,:)+v_coord(181,:));
        v_coord(233,:) = 0.25*(v_coord(16,:)+v_coord(237,:)+v_coord(173,:)+v_coord(53,:));
        v_coord(241,:) = 0.25*(v_coord(237,:)+v_coord(31,:)+v_coord(181,:)+v_coord(61,:));
        
        v_coord(83,:) = 0.5*(v_coord(38,:)+v_coord(113,:));
        v_coord(87,:) = 0.5*(v_coord(42,:)+v_coord(117,:));
        v_coord(91,:) = 0.5*(v_coord(46,:)+v_coord(121,:));  
        
        v_coord(111,:) = 0.5*(v_coord(8,:)+v_coord(113,:));
        v_coord(115,:) = 0.5*(v_coord(113,:)+v_coord(117,:));
        v_coord(119,:) = 0.5*(v_coord(117,:)+v_coord(121,:));
        v_coord(123,:) = 0.5*(v_coord(121,:)+v_coord(23,:));
        
        v_coord(143,:) = 0.5*(v_coord(113,:)+v_coord(173,:));
        v_coord(147,:) = 0.5*(v_coord(117,:)+v_coord(177,:));
        v_coord(151,:) = 0.5*(v_coord(121,:)+v_coord(181,:));  
        
        v_coord(171,:) = 0.5*(v_coord(12,:)+v_coord(173,:));
        v_coord(175,:) = 0.5*(v_coord(173,:)+v_coord(177,:));
        v_coord(179,:) = 0.5*(v_coord(177,:)+v_coord(181,:));
        v_coord(183,:) = 0.5*(v_coord(181,:)+v_coord(27,:));
        
        v_coord(203,:) = 0.5*(v_coord(173,:)+v_coord(233,:));
        v_coord(207,:) = 0.5*(v_coord(177,:)+v_coord(237,:));
        v_coord(211,:) = 0.5*(v_coord(181,:)+v_coord(241,:));  
        
        v_coord(231,:) = 0.5*(v_coord(16,:)+v_coord(233,:));       
        v_coord(235,:) = 0.5*(v_coord(233,:)+v_coord(237,:));
        v_coord(239,:) = 0.5*(v_coord(237,:)+v_coord(241,:));
        v_coord(243,:) = 0.5*(v_coord(241,:)+v_coord(31,:));
        
        v_coord(263,:) = 0.5*(v_coord(233,:)+v_coord(53,:));
        v_coord(267,:) = 0.5*(v_coord(237,:)+v_coord(57,:));
        v_coord(271,:) = 0.5*(v_coord(241,:)+v_coord(61,:));  
        
        v_coord(81,:) = 0.25*(v_coord(6,:)+v_coord(83,:)+v_coord(36,:)+v_coord(111,:));
        v_coord(85,:) = 0.25*(v_coord(83,:)+v_coord(87,:)+v_coord(40,:)+v_coord(115,:));
        v_coord(89,:) = 0.25*(v_coord(87,:)+v_coord(91,:)+v_coord(119,:)+v_coord(44,:));
        v_coord(93,:) = 0.25*(v_coord(91,:)+v_coord(21,:)+v_coord(48,:)+v_coord(123,:));

        v_coord(141,:) = 0.25*(v_coord(10,:)+v_coord(143,:)+v_coord(111,:)+v_coord(171,:));
        v_coord(145,:) = 0.25*(v_coord(143,:)+v_coord(147,:)+v_coord(115,:)+v_coord(175,:));
        v_coord(149,:) = 0.25*(v_coord(147,:)+v_coord(151,:)+v_coord(119,:)+v_coord(179,:));
        v_coord(153,:) = 0.25*(v_coord(151,:)+v_coord(25,:)+v_coord(123,:)+v_coord(183,:));
        
        v_coord(201,:) = 0.25*(v_coord(14,:)+v_coord(203,:)+v_coord(171,:)+v_coord(231,:));
        v_coord(205,:) = 0.25*(v_coord(203,:)+v_coord(207,:)+v_coord(175,:)+v_coord(235,:));
        v_coord(209,:) = 0.25*(v_coord(207,:)+v_coord(211,:)+v_coord(179,:)+v_coord(239,:));
        v_coord(213,:) = 0.25*(v_coord(211,:)+v_coord(29,:)+v_coord(183,:)+v_coord(243,:));
        
        v_coord(261,:) = 0.25*(v_coord(18,:)+v_coord(263,:)+v_coord(231,:)+v_coord(51,:));
        v_coord(265,:) = 0.25*(v_coord(263,:)+v_coord(267,:)+v_coord(235,:)+v_coord(55,:));
        v_coord(269,:) = 0.25*(v_coord(267,:)+v_coord(271,:)+v_coord(239,:)+v_coord(59,:));
        v_coord(273,:) = 0.25*(v_coord(271,:)+v_coord(33,:)+v_coord(243,:)+v_coord(63,:));
        
        for i = 0:6
            v_coord(66+2*i,:) = 0.5*v_coord(36+2*i,:)+0.5*v_coord(81+2*i,:);
            v_coord(96+2*i,:) = 0.5*v_coord(81+2*i,:)+0.5*v_coord(111+2*i,:);
            v_coord(126+2*i,:) = 0.5*v_coord(111+2*i,:)+0.5*v_coord(141+2*i,:);
            v_coord(156+2*i,:) = 0.5*v_coord(141+2*i,:)+0.5*v_coord(171+2*i,:);
            v_coord(186+2*i,:) = 0.5*v_coord(171+2*i,:)+0.5*v_coord(201+2*i,:);
            v_coord(216+2*i,:) = 0.5*v_coord(201+2*i,:)+0.5*v_coord(231+2*i,:);
            v_coord(246+2*i,:) = 0.5*v_coord(231+2*i,:)+0.5*v_coord(261+2*i,:);
            v_coord(276+2*i,:) = 0.5*v_coord(261+2*i,:)+0.5*v_coord(51+2*i,:);
        end
        
        v_coord(80,:) = 0.5*v_coord(6,:)+0.5*v_coord(81,:);
        v_coord(110,:) = 0.5*v_coord(8,:)+0.5*v_coord(111,:);
        v_coord(140,:) = 0.5*v_coord(10,:)+0.5*v_coord(141,:);
        v_coord(170,:) = 0.5*v_coord(12,:)+0.5*v_coord(171,:);
        v_coord(200,:) = 0.5*v_coord(14,:)+0.5*v_coord(201,:);
        v_coord(230,:) = 0.5*v_coord(16,:)+0.5*v_coord(231,:);
        v_coord(260,:) = 0.5*v_coord(18,:)+0.5*v_coord(261,:);
        for i = 0:5
            v_coord(82+2*i,:) = 0.5*v_coord(81+2*i,:)+0.5*v_coord(83+2*i,:);
            v_coord(112+2*i,:) = 0.5*v_coord(111+2*i,:)+0.5*v_coord(113+2*i,:);
            v_coord(142+2*i,:) = 0.5*v_coord(141+2*i,:)+0.5*v_coord(143+2*i,:);
            v_coord(172+2*i,:) = 0.5*v_coord(171+2*i,:)+0.5*v_coord(173+2*i,:);
            v_coord(202+2*i,:) = 0.5*v_coord(201+2*i,:)+0.5*v_coord(203+2*i,:);
            v_coord(232+2*i,:) = 0.5*v_coord(231+2*i,:)+0.5*v_coord(233+2*i,:);
            v_coord(262+2*i,:) = 0.5*v_coord(261+2*i,:)+0.5*v_coord(263+2*i,:);
        end
        v_coord(94,:) = 0.5*v_coord(93,:)+0.5*v_coord(21,:);
        v_coord(124,:) = 0.5*v_coord(123,:)+0.5*v_coord(23,:);
        v_coord(154,:) = 0.5*v_coord(153,:)+0.5*v_coord(25,:);
        v_coord(184,:) = 0.5*v_coord(183,:)+0.5*v_coord(27,:);
        v_coord(214,:) = 0.5*v_coord(213,:)+0.5*v_coord(29,:);
        v_coord(244,:) = 0.5*v_coord(243,:)+0.5*v_coord(31,:);
        v_coord(274,:) = 0.5*v_coord(273,:)+0.5*v_coord(33,:);

        v_coord(65,:) = 0.25*(v_coord(5,:)+v_coord(66,:)+v_coord(35,:)+v_coord(80,:));
        v_coord(95,:) = 0.25*(v_coord(7,:)+v_coord(96,:)+v_coord(80,:)+v_coord(110,:));
        v_coord(125,:) = 0.25*(v_coord(9,:)+v_coord(126,:)+v_coord(110,:)+v_coord(140,:));
        v_coord(155,:) = 0.25*(v_coord(11,:)+v_coord(156,:)+v_coord(140,:)+v_coord(170,:));
        v_coord(185,:) = 0.25*(v_coord(13,:)+v_coord(186,:)+v_coord(170,:)+v_coord(200,:));
        v_coord(215,:) = 0.25*(v_coord(15,:)+v_coord(216,:)+v_coord(200,:)+v_coord(230,:));
        v_coord(245,:) = 0.25*(v_coord(17,:)+v_coord(246,:)+v_coord(230,:)+v_coord(260,:));
        v_coord(275,:) = 0.25*(v_coord(19,:)+v_coord(276,:)+v_coord(260,:)+v_coord(50,:));
        for i = 0:5
            v_coord(67+2*i,:) = 0.25*(v_coord(66+2*i,:)+v_coord(68+2*i,:)+v_coord(37+2*i,:)+v_coord(82+2*i,:));
            v_coord(97+2*i,:) = 0.25*(v_coord(96+2*i,:)+v_coord(98+2*i,:)+v_coord(82+2*i,:)+v_coord(112+2*i,:));
            v_coord(127+2*i,:) = 0.25*(v_coord(126+2*i,:)+v_coord(128+2*i,:)+v_coord(112+2*i,:)+v_coord(142+2*i,:));
            v_coord(157+2*i,:) = 0.25*(v_coord(156+2*i,:)+v_coord(158+2*i,:)+v_coord(142+2*i,:)+v_coord(172+2*i,:));
            v_coord(187+2*i,:) = 0.25*(v_coord(186+2*i,:)+v_coord(188+2*i,:)+v_coord(172+2*i,:)+v_coord(202+2*i,:));
            v_coord(217+2*i,:) = 0.25*(v_coord(216+2*i,:)+v_coord(218+2*i,:)+v_coord(202+2*i,:)+v_coord(232+2*i,:));
            v_coord(247+2*i,:) = 0.25*(v_coord(246+2*i,:)+v_coord(248+2*i,:)+v_coord(232+2*i,:)+v_coord(262+2*i,:));
            v_coord(277+2*i,:) = 0.25*(v_coord(276+2*i,:)+v_coord(278+2*i,:)+v_coord(262+2*i,:)+v_coord(52+2*i,:));
        end
        v_coord(79,:) = 0.25*(v_coord(78,:)+v_coord(20,:)+v_coord(49,:)+v_coord(94,:));
        v_coord(109,:) = 0.25*(v_coord(108,:)+v_coord(22,:)+v_coord(94,:)+v_coord(124,:));
        v_coord(139,:) = 0.25*(v_coord(138,:)+v_coord(24,:)+v_coord(124,:)+v_coord(154,:));
        v_coord(169,:) = 0.25*(v_coord(168,:)+v_coord(26,:)+v_coord(154,:)+v_coord(184,:));
        v_coord(199,:) = 0.25*(v_coord(198,:)+v_coord(28,:)+v_coord(184,:)+v_coord(214,:));
        v_coord(229,:) = 0.25*(v_coord(228,:)+v_coord(30,:)+v_coord(214,:)+v_coord(244,:));
        v_coord(259,:) = 0.25*(v_coord(258,:)+v_coord(32,:)+v_coord(244,:)+v_coord(274,:));
        v_coord(289,:) = 0.25*(v_coord(288,:)+v_coord(34,:)+v_coord(274,:)+v_coord(64,:));
    end
    
end
  
v_coord
plot(v_coord(:,1),v_coord(:,2),'bo')
end