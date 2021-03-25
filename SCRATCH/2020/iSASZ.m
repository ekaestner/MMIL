ttt = summer(16);
[ rgb('greenish grey') ; ttt([10 14],:) ]

[ rgb2hsv( rgb('greenish grey') ) ;
  rgb2hsv( ttt([10],:) ) ;
  rgb2hsv( ttt([14],:)) ]

[ hsv2rgb( [ .29 .02 .50 ] )
  hsv2rgb( [ .29 .20 .40 ] ) ;
  hsv2rgb( [ .27 .40 .50 ] ) ; 
  hsv2rgb( [ .23 .65 .70 ] ) ;
  hsv2rgb( [ .17 .95 .90 ] ) ]


lhs_ovr_ele_loc( string_find( lhs_ovr_ele_loc(:,1) , 'NY007' ) ,:) = [];

%%
'dark blue' 'blue' 'bright blue'

[ rgb2hsv( rgb('dark blue') ) ;
  rgb2hsv( rgb('blue') ) ;
  rgb2hsv( rgb('bright blue') ) ]

%%
'dark red' 'red' 'bright red'
[ rgb2hsv( rgb('dark red') ) ;
  rgb2hsv( rgb('red') ) ;
  rgb2hsv( rgb('bright red') ) ]

%%
'dark purple' 'purple' 'bright purple'
[ rgb2hsv( rgb('dark purple') ) ;
  rgb2hsv( rgb('purple') ) ;
  rgb2hsv( rgb('bright purple') ) ]

%%
tot_cnt.visual( cellfun( @isempty , tot_cnt.visual ) ) = {nan}
tot_cnt.auditory( cellfun( @isempty , tot_cnt.auditory ) ) = {nan}
tot_cnt.tot( cellfun( @isempty , tot_cnt.tot ) ) = {nan}

cell2mat(tot_cnt.visual(2:9,2:9)) ./ cell2mat(tot_cnt.tot(2:9,2:9))

hld = tot_cnt.visual;
hld(2:9,2:9) = num2cell( round((cell2mat(tot_cnt.visual(2:9,2:9)) ./ cell2mat(tot_cnt.tot(2:9,2:9))) * 100) );
hld( [1 3 6 7 8] , [1 3 6 7 8 ] )

hld = tot_cnt.auditory;
hld(2:9,2:9) = num2cell( round((cell2mat(tot_cnt.auditory(2:9,2:9)) ./ cell2mat(tot_cnt.tot(2:9,2:9))) * 100) );
hld( [1 3 6 7 8] , [1 3 6 7 8 ] )