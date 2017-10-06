	!	---SCI_EX_02.F90
	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW LINE_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(num_RANGES=3,num_AXIS=2,num_perange=4,num_Labels=4, &
		num_legends=3 )
	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit

	RECORD /GraphSettings/ lineGraph
	RECORD /DataSettings/ lineDataSets(num_RANGES)  ! 3 data sets (ranges)
	RECORD /AxisSettings/ lineAxes(num_AXIS)        ! 2 axes

	REAL*4      lineData(num_RANGES,num_perange)    ! 3 data sets each with 4 (y)
	CHARACTER*20 lineLabels(num_Labels)             ! the x labels for the line  graph
	CHARACTER*20 lineDataLegends(num_legends)       ! 3 data range legends

	INTEGER   retcode
	INTEGER   setLegends
	RECORD /windowconfig/ wc
	!	TO SET THE DATA SETS FOR YOUR GRAPH.
	DATA lineData /                       &
		300.0, 350.4, 402.5, 380.1,    &
		450.6, 395.7, 410.2, 450.9,    &
		500.2, 550.1, 520.4, 590.8/

	DATA lineLabels / 'FIRST','SECOND','THIRD','FOURTH' /

	DATA lineDataLegends / '1995','1996','1997'/

	i_quit = setexitqq(QWIN$EXITPERSIST)

	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!	CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTLINE,lineGraph)
	lineGraph.setGraphMode=.FALSE.
	lineGraph.x2=wc.numxpixels-1
	lineGraph.y2=wc.numypixels-1
	lineGraph.graphbgcolor=$ciwhite
	lineGraph.title='General Production'
	lineGraph.titlecolor=$ciblack
	lineGraph.title2='(Million Yuan)'
	lineGraph.title2color=$ciBlack
	lineGraph.graphcolor=$ciBlack

	!	THEN, CALL 	GetLabelMultiDataDefaults() TO SET LABEL DATA.
	retcode=GetLabelMultiDataDefaults(lineGraph,num_Labels,lineLabels, &
		lineData,num_legends,lineDataSets)
	DO setLegends=1,num_legends
		lineDataSets(setLegends).title=lineDataLegends(setLegends)
	END DO
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(lineGraph,num_labels,lineDataSets, &
		$ATX,$AFLINEAR,lineAxes(1))
	lineAxes(1).title='Quater of the Year'

	retcode=GetAxisMultiDefaults(lineGraph,num_Labels,lineDataSets, &
		$ATY,$AFLINEAR,lineAxes(2))
	lineAxes(2).title='General Production'
	lineAxes(2).tickType=$TTOUTSIDE
	lineAxes(2).numDigits=1
	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(lineGraph,num_axis,lineAxes,num_legends)
	!   FINALLY, CALL PlotLabelMultiData() TO WRITE LABELS.
	retcode=PlotLabelMultiData(lineGraph,lineLabels,lineData, &
		num_ranges,lineDataSets,lineAxes(1),lineAxes(2))
	END
