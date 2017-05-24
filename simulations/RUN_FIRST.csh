
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

	while( $s < 101 )
		foreach scen ( 1 2 3 4 5 )
			# par1: $s seed number
			# par2: 200 the sample size n
		  # par3: 6 or 100, number of variables
		  # par4: 1 the noise variance
		  # par5: $scen, the scenario number

			qsub -q "shojaie*"  -cwd  -o /dev/null -e /dev/null  simulations.csh $s 200 6 1 $scen
		end
		@ s++
	end




