
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

	while( $s < 31 )
		# par1: $s seed number
		# par2: nvar 100; number of variabls to use.

		qsub -q "shojaie*"  -cwd  -o /dev/null -e /dev/null  simulations.csh $s 100
		
		@ s++
	end




