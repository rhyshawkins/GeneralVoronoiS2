
TARGETS = synthetic_template.txt \
	synthetic_obs_constant.txt \
	results_constant/ch.dat \
	results_constant/mean.txt 

all : $(TARGETS)

synthetic_template.txt :
	python2 ../../scripts/mksyntheticpointstemplate.py -o $@ -n 100

synthetic_obs_constant.txt : synthetic_template.txt
	../mksyntheticregression -i $< -o $@ \
	-m Constant \
	-n 0.5 \
	-O synthetic_obs_constant.true \
	-s 10.0 -f 0.0

prior_values.txt :
	echo "prior Uniform" > $@
	echo "0.0 10.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@

prior_position.txt :
	echo "sphericalprior UniformSpherical" > $@
	echo "sphericalproposal VonMisesSpherical" >> $@
	echo "1.0e2" >> $@

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


