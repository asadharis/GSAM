
#!/bin/csh

#Run a loop to do our analysis for each seed
set s = 1

while( $s < 501 )
		foreach scen ( 1 2 3 4 5 )
		  foreach n ( 100 200 300 400 500 600 700 800 )

			  # par1: $s seed number
			  # par2: 200 the sample size n
		    # par3: 6 or 100, number of variables
		    # par4: 1 the noise variance
		    # par5: $scen, the scenario number

			  qsub -q "shojaie*"  -cwd  -o /dev/null -e /dev/null  simulations.csh $s $n 6 1 $scen
			end
		end
	@ s++
end




