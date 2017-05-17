
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

	while( $s < 101 )
		foreach scen ( 1 2 3 4 5 ) 
	
			qsub -q "shojaie*"  -cwd  -o /dev/null -e /dev/null  simulations.csh $s 200 6 10 $scen
		end
		@ s++
	end




