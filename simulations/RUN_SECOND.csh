
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

	while( $s < 501  )
		foreach scen ( 1 2 3 4 5 )
		  # par1: $s seed number
		  # par2: 6 or 100, number of variables
		  # par3: 1 the noise variance
		  # par4: $scen, the scenario number

			qsub -q "shojaie*"  -cwd  -o /dev/null -e /dev/null  simulations2.csh $s 6 1 $scen
		end
		@ s++
	end




