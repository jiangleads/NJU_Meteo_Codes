	!	---SCI_EX_04.F90
	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW XY_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(num_RANGES=1,num_AXIS=4,num_perange=400,num_Labels=4)
	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit
	RECORD /GraphSettings/ xyGraph
	RECORD /DataSettings/ xyDataSets(num_RANGES)    ! number of data sets (ranges)
	RECORD /AxisSettings/ xyAxes(num_AXIS)        ! number of  axes

	REAL*4      xyData(2,num_perange,num_RANGES) !  data sets

	INTEGER     retcode
	RECORD /windowconfig/ wc
	!	TO SET THE DATA SETS FOR YOUR GRAPH.

	DO index=1,num_perange
		theta=$PI*(index-1.0D0)/(num_perange-1.0D0)
		xyData(1,index,1)=2*sin(5*theta)
		xyData(2,index,1)=theta
	END DO

	i_quit = setexitqq(QWIN$EXITPERSIST)
	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!   CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTXY,xyGraph)
	xyGraph.setGraphMode=.FALSE.
	xyGraph.x2=wc.numxpixels-1
	xyGraph.y2=wc.numypixels-1
	xyGraph.title='SINE GRAPH'
	xyGraph.titlecolor=$ciblack
	xyGraph.graphbgcolor=$ciwhite
	xyGraph.graphcolor=$ciblack
	!	THEN, CALL 	GetLabelMultiDataDefaults() TO SET LABEL DATA.
	retcode=GetMultiDataDefaults(xyGraph,num_perange,xyData, &
		num_RANGES,xyDataSets)

	xyDataSets(1).linecolor=$cired
	xyDataSets(1).markertype=0
	xyDataSets(1).linetype=$LTthicksolid
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(xyGraph,num_ranges,xyDataSets,  &
		$ATX,$AFLINEAR,xyAxes(1))
	xyAxes(1).title='THETA'

	retcode=GetAxisMultiDefaults(xyGraph,num_ranges,xyDataSets,  &
		$ATY,$AFLINEAR,xyAxes(2))

	retcode=GetAxisMultiDefaults(xyGraph,num_RANGES,xyDataSets,  &
		$ATX,$AFLINEAR,xyAxes(3))
	xyAxes(3).lowVal=xyAxes(3).lowVal/2.54
	xyAxes(3).highVal=xyAxes(3).highVal/2.54
	xyAxes(3).increment=xyAxes(3).increment/2.54
	xyAxes(3).tickColor=$CIWHITE
	xyAxes(3).minorTickColor=$CIGRAY

	retcode=GetAxisMultiDefaults(xyGraph,num_RANGES,xyDataSets,  &
		$ATY,$AFLINEAR,xyAxes(4))

	xyAxes(4).lowVal=xyAxes(4).lowVal/2.54
	xyAxes(4).highVal=xyAxes(4).highVal/2.54
	xyAxes(4).increment=xyAxes(4).increment/2.54
	xyAxes(4).tickColor=$CIWHITE
	xyAxes(4).minorTickColor=$CIGRAY

	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(xyGraph,num_AXIS,xyAxes,num_RANGES)
	!   FINALLY, CALL PlotLabelMultiData() TO WRITE LABELS.
	retcode=PlotMultiData(xyGraph,xyData,num_RANGES, &
		xyDataSets,xyAxes(1),xyAxes(2))

	END
