
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

	while( $s < 2 )
		foreach scen ( 1 )#2 3 4 5 ) 
	
			qsub -q shojaie-normal.q  -cwd  -o /dev/null -e /dev/null  simulations2.csh $s 6 10 $scen
		end
		@ s++
	end




