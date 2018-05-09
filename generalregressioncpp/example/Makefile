
TARGETS = synthetic_template.txt \
	synthetic_obs_constant.txt \
	results_constant/ch.dat \
	results_constant/mean.txt 

all : $(TARGETS)

synthetic_template.txt :
	python2 ../../scripts/mksyntheticpointstemplate.py -o $@ -n 1000

synthetic_obs_constant.txt : synthetic_template.txt
	../mksyntheticregression -i $< -o $@ \
	-m Constant \
	-n 0.5 \
	-O synthetic_obs_constant.true \
	-s 10.0 -f 0.0

synthetic_obs_checkerboard.txt : synthetic_template.txt
	../mksyntheticregression -i $< -o $@ \
	-m CheckerBoard \
	-n 0.5 \
	-O synthetic_obs_constant.true \
	-I synthetic_obs_checkerboard.img \
	-s 10.0 -f 0.0

prior_values.txt :
	echo "prior Uniform" > $@
	echo "0.0 10.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.1" >> $@

prior_position.txt :
	echo "sphericalprior UniformSpherical" > $@
	echo "sphericalproposal VonMisesSpherical" >> $@
	echo "2.0e2" >> $@

prior_hierarchical.txt :
	echo "prior Uniform" > $@
	echo "0.0 10.0" >> $@
	echo "proposal LogGaussian" >> $@
	echo "0.1" >> $@


results_constant/ch.dat : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt
	mkdir -p results_constant
	../regressionvoronois2 -i synthetic_obs_constant.txt -o results_constant/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 10 \
	-t 50000

results_constant/mean.txt : results_constant/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 10000 \
	-z 0 -Z 10 \
	-m results_constant/mean.txt \
	-T results_constant/stddev.txt \
	-W 180 -H 90

results_checkerboard/ch.dat : prior_values.txt \
	prior_position.txt \
	prior_hierarchical.txt \
	synthetic_obs_checkerboard.txt
	mkdir -p results_checkerboard
	../regressionvoronois2 -i synthetic_obs_checkerboard.txt -o results_checkerboard/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-H prior_hierarchical.txt \
	-T 100 \
	-t 1000000 \
	-v 10000

results_checkerboard/mean.txt : results_checkerboard/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 100 -s 500000 \
	-z 0 -Z 10 \
	-m results_checkerboard/mean.txt \
	-T results_checkerboard/stddev.txt \
	-W 90 -H 45


