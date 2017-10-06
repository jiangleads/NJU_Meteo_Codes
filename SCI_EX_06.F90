	!	---SCI_EX_06.F90
	!	---THIS PROGRAM USES SCIGRAPH PACKAGE TO DRAW POLAR_GRAPH.
	!
	USE MSFLIB
	USE SCIGRAPH
	PARAMETER(num_RANGES=2,num_AXIS=2,num_p=200,NUM_LEGENDS=2)

	!	---TO DIMENSION VARIABLES.
	INTEGER i_quit

	RECORD /GraphSettings/ polarGraph
	RECORD /DataSettings/ polarDataSets(num_RANGES)
	RECORD /AxisSettings/ polarAxes(num_AXIS)


	REAL*4      polarData(2,NUM_P,NUM_RANGES)! 2 data set each with NUM_P (r,th)
	CHARACTER*20 polarDataLegends(NUM_LEGENDS)! 2 data range legends

	INTEGER     retcode
	INTEGER     index
	REAL*4      theta

	RECORD /windowconfig/ wc

	DATA polarDataLegends / 'r=2Sin(5@)','r=1-2cos(@)'/
	i_quit = setexitqq(QWIN$EXITPERSIST)
	if( .not. GetWindowConfig(wc) ) stop 'Window Not Open'

	! set up data range 1
	DO index=1,num_p
		theta=$PI*(index-1.0D0)/(num_p-1.0D0)
		polarData(1,index,1)=2*sin(5*theta)
		polarData(2,index,1)=theta
	END DO

	! set up data range 2
	DO index=1,num_p
		theta=2*$PI*(index-1.0)/(num_p-1.0)
		polarData(1,index,2)=1-2*cos(5*theta)
		polarData(2,index,2)=theta
	END DO
	!	WHEN YOU PLOT GRAPH WITH SCIGRAPH PACKAGE,
	!   CALL GetGraphDefaults() first.
	retcode=GetGraphDefaults($GTPOLAR,polarGraph)
	polarGraph.setGraphMode=.FALSE.
	polarGraph.x2=wc.numxpixels-1
	polarGraph.y2=wc.numypixels-1
	polarGraph.title='极坐标图示例'
	polarGraph.graphbgcolor=$ciwhite
	!	THEN, CALL 	GetMultiDataDefaults() TO SET up DATA.
	retcode=GetMultiDataDefaults(polarGraph,num_p, &
		polarData,NUM_RANGES,polarDataSets)
	polarDataSets(1).title=polarDataLegends(1)
	polarDataSets(1).lineColor=$CIblack
	polarDataSets(1).markerType=$MKNONE

	polarDataSets(2).title=polarDataLegends(2)
	polarDataSets(2).lineColor=$CIlightblue
	polarDataSets(2).markerType=$MKNONE
	!	AND THEN, CALL GetAxisMultiDefaults() TO SET AXIS DATA.
	retcode=GetAxisMultiDefaults(polarGraph,NUM_RANGES,polarDataSets, &
		$ATR,$AFLINEAR,polarAxes(1))

	retcode=GetAxisMultiDefaults(polarGraph,NUM_RANGES,polarDataSets, &
		$ATTHETA,$AFLINEAR,polarAxes(2))
	polarAxes(2).gridStyle=$GSBOTH
	polarAxes(2).gridLineType=$LTDOT
	polarAxes(2).gridColor=$CIblack
	!	FOURTHLY, CALL 	PlotGraph() TO DRAW THE GRAPH.
	retcode=PlotGraph(polarGraph,NUM_AXIS,polarAxes,NUM_RANGES)
	!   FINALLY, CALL PlotMultiData() TO COMPLETE THE GRAPH.
	retcode=PlotMultiData(polarGraph,polarData,         &
		NUM_RANGES,polarDataSets,polarAxes(1),polarAxes(2))

	END
