
TARGETS = synthetic_template.txt \
	synthetic_obs_constant.txt \
	results_constant/ch.dat-000 \
	results_constant/ch.dat-001 \
	results_constant/mean.txt 

all : $(TARGETS)

#
# Generate a random set of paths using the provided script
#
synthetic_template.txt :
	python2 ../../scripts/mksyntheticpathtemplate.py -o $@ -n 20

#
# Create observations using a constant model. The mksynthetictomography program will generate a
# normalized model with a constant 0.5 which will become 3.0 with the
# scale and offset parameters -s and -f, ie 0.5 * 3.0 + 1.5 = 3.0
# The template generates surface paths at Earth radii with corresponding
# times of around 2000s, so noise level of 10.0s corresponds to approximately
# a 0.5 percent error.
#
synthetic_obs_constant.txt : synthetic_template.txt
	../mksynthetictomography -i $< -o $@ \
	-m Constant \
	-n 10.0 \
	-O synthetic_obs_constant.true \
	-s 3.0 -f 1.5

#
# In this synthetic example, with the scale and offset parameters
# above the maximum range is 1.5 to 4.5 km/s so we've used a uniform
# prior about this range. The proposal is Gaussian with a small
# standard deviation of 0.005. This last value should be tuned to
# ensure an adequate acceptance rate of the Value proposals, ie
# between 20 and 50 percent.
#
prior_values.txt :
	echo "prior Uniform" > $@
	echo "1.5 4.5" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.005" >> $@

prior_position.txt :
	echo "sphericalprior UniformSpherical" > $@
	echo "sphericalproposal VonMisesSpherical" >> $@
	echo "1.0e2" >> $@

#
# For synthetic studies where we know the true noise we don't need
# to use hierarchical sampling. However adding it in gives confidence
# that we're able to recover model + noise in an inversion. This
# specifies the prior and proposal for the hierachical scaling
# parameter. This is a generally good prior/proposal combination
# where the use of a truncated Jeffrey's prior means we are
# scale invariant and the LogGaussian proposal perturbs the
# parameter with a Gaussian in log-space meaning the perturbation
# is of an appropriate scale.
#
prior_hierarchical.txt :
	echo "prior Jeffreys" > $@
	echo "1.0e-9 1.0e1" >> $@
	echo "proposal LogGaussian" >> $@
	echo "1.0e-2" >> $@

#
# Running an inversion produces a chain history file (ch.dat)
#
results_constant/ch.dat-000 : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt
	mkdir -p results_constant
	mpirun -np 1 ../tomographyvoronois2pt -i synthetic_obs_constant.txt -o results_constant/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000 -v 1000 \
	-O 0

results_constant/ch.dat-001 : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt
	mkdir -p results_constant
	mpirun -np 1 ../tomographyvoronois2pt -i synthetic_obs_constant.txt -o results_constant/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000 -v 1000 \
	-O 1

#
# To output images, we use post processing
#
results_constant/mean.txt : results_constant/ch.dat-000 results_constant/ch.dat-001
	mpirun -np 2 ../../postS2Voronoi_mean_mpi -i results_constant/ch.dat -o $@ \
	-t 10 -s 25000 \
	-z 1.5 -Z 4.5 \
	-m results_constant/median.txt \
	-T results_constant/stddev.txt \
	-W 180 -H 90

results_constant_cont/ch.dat-000 : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt \
	results_constant/ch.dat-000
	mkdir -p results_constant_cont
	mpirun -np 1 ../tomographyvoronois2pt -i synthetic_obs_constant.txt -o results_constant_cont/ \
	-I results_constant/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000 -v 1000 \
	-O 0 

results_constant_cont/ch.dat-001 : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt \
	results_constant/ch.dat-001
	mkdir -p results_constant_cont
	mpirun -np 1 ../tomographyvoronois2pt -i synthetic_obs_constant.txt -o results_constant_cont/ \
	-I results_constant/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000 -v 1000 \
	-O 1

#
# To output history of the likelihood/hierarchical/no. cells(k) to monitor convergence
# Above inversion wasn't hierarchical. 
#
results_constant/like.txt : results_constant/ch.dat
	../../postS2Voronoi_likelihood -i $< -o $@ \
	-K results_constant/khistory.txt \
	-t 10 



#
# Same paths but observations with a East/West hemispherical model
#
synthetic_obs_ew.txt : synthetic_template.txt
	../mksynthetictomography -i $< -o $@ \
	-m EastWest \
	-n 10.0 \
	-O synthetic_obs_ew.true \
	-s 3.0 -f 1.5

#
# Same inversion but observations with a East/West hemispherical model
#
results_eastwest/ch.dat : prior_values.txt \
	prior_position.txt \
	synthetic_obs_ew.txt
	mkdir -p results_eastwest
	../tomographyvoronois2 -i synthetic_obs_ew.txt -o results_eastwest/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000 -v 100

results_eastwest/mean.txt : results_eastwest/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 25000 \
	-z 1.5 -Z 4.5 \
	-m results_eastwest/median.txt \
	-T results_eastwest/stddev.txt \
	-W 180 -H 90

#
# Inversion using hierarchical error sampling (note the extra -H prior_hierarchical.txt
# option.
#
results_eastwest_hierarchical/ch.dat : prior_values.txt \
	prior_position.txt \
	prior_hierarchical.txt \
	synthetic_obs_ew.txt
	mkdir -p results_eastwest_hierarchical
	../tomographyvoronois2 -i synthetic_obs_ew.txt -o results_eastwest_hierarchical/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-H prior_hierarchical.txt \
	-T 10 \
	-t 50000 -v 1000

#
# Image is post processed in the same way
#
results_eastwest_hierarchical/mean.txt : results_eastwest_hierarchical/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 25000 \
	-z 1.5 -Z 4.5 \
	-m results_eastwest_hierarchical/median.txt \
	-T results_eastwest_hierarchical/stddev.txt \
	-W 180 -H 90

#
# Extra hierarchical information 
#
results_eastwest_hierarchical/like.txt : results_eastwest_hierarchical/ch.dat
	../../postS2Voronoi_likelihood -i $< -o $@ \
	-H results_eastwest_hierarchical/hierarchical.txt \
	-K results_eastwest_hierarchical/khistory.txt \
	-t 10 
