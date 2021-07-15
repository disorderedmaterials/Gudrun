	Integer*4 ngroupt			!defines no. of groups formed
	integer*4 nq			!number of Q values
	integer*4 nweighterr		!1 means to use error bars in merge, else 0
	integer*4 nqfirst		!gives the first element of the final merged data that is to be output
	real*4 aggweights(mq,mgroup)	!aggregated weights array
	real*4 aggsweights(mq,mgroup)	!aggregated sample weights array
	real*4 aggeweights(mq,mgroup)	!aggregated error weights array
	real*4 aggdweights(mq,mgroup)	!aggregated deviation weights array
	real*4 Qbound(mq)		!q boundary values
	real*4 qvalue(mq)		!q values
	common/merge_arrays_i/ngroupt,nq,nweighterr,nqfirst
	common/merge_arrays_r/aggweights,aggsweights,aggeweights
     $,aggdweights,qbound,qvalue

