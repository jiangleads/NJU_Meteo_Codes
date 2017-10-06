	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW BAR_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(num_RANGES=3,num_AXIS=2,num_perange=4,num_Labels=4, &
		num_legends=3 )
	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit

	RECORD /GraphSettings/ barGraph
	RECORD /DataSettings/ barDataSets(num_RANGES)    ! number of data sets (ranges)
	RECORD /AxisSettings/ barAxes(num_AXIS)     ! number of axes

	REAL*4      barData(num_RANGES,num_perange)
	CHARACTER*20 barLabels(num_Labels)        ! the x labels for the bar graph
	CHARACTER*20 barDataLegends(num_legends)  ! number of data range legends

	INTEGER     retcode
	INTEGER     setLegends
	RECORD /windowconfig/ wc
	!	TO SET THE DATA SETS FOR YOUR GRAPH.
	DATA barData /                     &
		300.0, 350.4, 402.5, 380.1,    &
		450.6, 395.7, 410.2, 450.9,    &
		500.2, 550.1, 520.4, 590.8/

	DATA barLabels / '第一季度','第二季度','第三季度','第四季度'/

	DATA barDataLegends / '1995年', '1996年', '1997年' /
	i_quit = setexitqq(QWIN$EXITPERSIST)
	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!   CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTBAR,barGraph)
	barGraph.setGraphMode=.FALSE.
	barGraph.x2=wc.numxpixels-1
	barGraph.y2=wc.numypixels-1
	barGraph.title='某工厂连续三年季平均产值'
	barGraph.title2='单位：百万元'
	barGraph.titlecolor=$ciblack
	barGraph.graphbgcolor=$ciwhite
	barGraph.graphcolor=$ciblack
	!	THEN, CALL 	GetLabelMultiDataDefaults() TO SET LABEL DATA.
	retcode=GetLabelMultiDataDefaults(barGraph,num_labels,barLabels, &
		barData,num_legends,barDataSets)
	DO setLegends=1,num_legends
		barDataSets(setLegends).title=barDataLegends(setLegends)
	END DO
	barDataSets(num_legends).barType=$BTHASHLEFT
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(barGraph,num_labels,barDataSets,  &
		$ATX,$AFLINEAR,barAxes(1))
	barAxes(1).title='季度'
	barAxes(1).titlecolor=$ciblack
	retcode=GetAxisMultiDefaults(barGraph,num_Labels,barDataSets,  &
		$ATY,$AFLINEAR,barAxes(2))
	barAxes(2).title='产值'
	barAxes(1).titlecolor=$ciblack
	!          barAxes(2).gridStyle=$GSBOTH
	!          barAxes(2).gridLineType=$LTDOT
	!          barAxes(2).gridColor=$CILIGHTBLUE
	barAxes(2).numDigits=1
	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(barGraph,num_axis,barAxes,num_legends)
	!   FINALLY, CALL PlotLabelMultiData() TO WRITE LABELS.
	retcode=PlotLabelMultiData(barGraph,barLabels,barData,  &
		num_ranges,barDataSets,barAxes(1),barAxes(2))
	END
