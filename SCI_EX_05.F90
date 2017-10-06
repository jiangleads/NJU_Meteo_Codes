	!	---SCI_EX_05.F90
	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW XY_ERROR_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(num_RANGES=2,num_AXIS=2,num_perange=7)
	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit

	RECORD /GraphSettings/ simplegraph
	RECORD /DataSettings/ simpledsets(num_RANGES)
	RECORD /AxisSettings/ simpleaxes(num_AXIS)

	REAL*4      simpledata(3,num_perange,num_RANGES)

	INTEGER     retcode
	RECORD /windowconfig/ wc
	!	TO SET THE DATA SETS FOR YOUR GRAPH.
	DATA simpledata / 1.0,2,0.15, 1.2,4,0.2, 1.3,6,0.3, 1.4,7,0,  &
		1.5,8,0.4, 1.9,9,0.125, 1.95,10,0.2,&
		1.0,4,0.3, 1.2,8,0.4, 1.3,12,0.6, 1.4,14,0.2,&
		1.5,12,0.8, 1.9,13,0.25, 1.95,15,0.4 /
	i_quit = setexitqq(QWIN$EXITPERSIST)
	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!   CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTXYWERRBAR,simplegraph)
	simplegraph.setGraphMode = .FALSE.
	simplegraph.x2 = wc.numxpixels-1
	simplegraph.y2 = wc.numypixels-1
	simplegraph.graphbgcolor =$ciwhite
	simplegraph.title = 'XY_ERROR GRAPH'
	simplegraph.titlecolor = $ciblack
	!	THEN, CALL 	GetMultiDataDefaults() TO SET  DATA.
	retcode=GetMultiDataDefaults(simplegraph,num_perange,&
		simpledata, num_RANGES,simpledsets)
	simpledsets.markerColor=$CIblack
	simpledsets.lineType=$LTsolid
	simpledsets.lineColor=$CIblue
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(simplegraph,num_RANGES,simpledsets,$ATX,  &
		$AFLINEAR,simpleaxes(1))
	simpleaxes(1).axisPos=$APTOP
	simpleaxes(1).titleColor=$CIblack
	simpleaxes(1).axisColor=$CIblack
	simpleaxes(1).tickColor=$CIblack
	simpleaxes(1).minorTickColor=$CIGREEN
	simpleaxes(1).tickType=$TTOUTSIDE

	retcode=GetAxisMultiDefaults(simplegraph,NUM_RANGES,simpledsets,$ATY,  &
		$AFLINEAR,simpleaxes(2))
	simpleaxes(2).titleColor=$CIblack
	simpleaxes(2).axisColor=$CIblack
	simpleaxes(2).tickColor=$CIblack
	simpleaxes(2).minorTickColor=$CIblack
	simpleaxes(2).tickType=$TTOUTSIDE
	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(simplegraph,NUM_AXIS,simpleaxes,NUM_PERANGE)
	!   FINALLY, CALL PLOTMultiData() TO COMPLETE THE GRAPH.
	retcode=PlotMultiData(simplegraph,simpledata,NUM_RANGES,simpledsets,  &
		simpleaxes(1),simpleaxes(2))

	END
