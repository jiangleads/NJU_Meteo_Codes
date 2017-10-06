	!	---SCI_EX_03.F90
	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW LINE_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(numRANGES=4,numAXIS=2,numperange=6,numLabels=6, &
		numlegends=4,numDim=2 )
	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit

	RECORD /GraphSettings/ lineGraph
	RECORD /DataSettings/ lineDataSets(numRANGES)  ! 4 data sets (ranges)
	RECORD /AxisSettings/ lineAxes(numAXIS)        ! 2 axes

	REAL*4      lineData(numperange,numRANGES, numDim)! 4 data sets each with 6 (y)
	CHARACTER*20 lineLabels(numLabels)               ! the x labels for the line graph
	CHARACTER*20 lineDataLegends(numlegends)          ! 4 data range legends

	INTEGER     retcode
	INTEGER     setLegends
	RECORD /windowconfig/ wc
	!	TO SET THE DATA SETS FOR YOUR GRAPH.

	DATA lineData /                                                    &
		316.5, 2, 317.8, 1, 320.1, 3, 318.3, 2.4, 314.2, 3, 315.1, 2.2,   &
		324.6, 0, 326.6, 9, 327.8, 3, 326.3, 1.2, 323.1, 5, 324.0, 2,     &
		337.8, 5, 339.9, 2, 341.2, 2, 339.3, .8,  335.7, 4, 336.7, 3,     &
		349.9, 1, 351.9, 4, 353.9, 5, 352.1, 1,   348.5, 2, 349.8, 1.1 /

	DATA lineLabels / 'January','March','May','July','September','November' /

	DATA lineDataLegends / '1960','1970','1980','1988'/
	i_quit = setexitqq(QWIN$EXITPERSIST)
	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!   CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTLINEWERRBAR,lineGraph)
	lineGraph.setGraphMode=.FALSE.
	lineGraph.x2=wc.numxpixels-1
	lineGraph.y2=wc.numypixels-1
	lineGraph.title='Atmospheric CO2'
	lineGraph.title2='Parts Per Million by Volume'
	lineGraph.graphbgcolor=$ciwhite
	lineGraph.graphcolor=$ciBlack
	!	THEN, CALL 	GetLabelMultiDataDefaults() TO SET LABEL DATA.
	retcode=GetLabelMultiDataDefaults(lineGraph,numperange,lineLabels,  &
		lineData,numranges,lineDataSets)
	DO setLegends=1,numlegends
		lineDataSets(setLegends).title=lineDataLegends(setLegends)
	END DO
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(lineGraph,numranges,lineDataSets,  &
		$ATX,$AFLINEAR,lineAxes(1))
	lineAxes(1).title='Month'

	retcode=GetAxisMultiDefaults(lineGraph,numranges,lineDataSets,  &
		$ATY,$AFLINEAR,lineAxes(2))
	lineAxes(2).title='Concentration of CO2'
	lineAxes(2).tickType=$TTOUTSIDE
	lineAxes(2).numDigits=1
	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(lineGraph,numAxis,lineAxes,numranges)
	!   FINALLY, CALL PlotLabelMultiData() TO WRITE LABELS.
	retcode=PlotLabelMultiData(lineGraph,lineLabels,lineData,  &
		numranges,lineDataSets,lineAxes(1),lineAxes(2))
	END
