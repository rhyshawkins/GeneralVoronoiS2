
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

results_constant_opt/ch.dat : prior_values.txt \
	prior_position.txt \
	synthetic_obs_constant.txt
	mkdir -p results_constant_opt
	../regressionvoronois2 -i synthetic_obs_constant.txt -o results_constant_opt/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 1000 \
	-t 1000 \
	-v 10 \
	-C 100 \
	-z 1000 \
	-Z 5.0

results_constant_opt/mean.txt : results_constant_opt/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 0 \
	-z 0 -Z 10 \
	-T results_constant_opt/stddev.txt \
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

results_checkerboard_opt/ch.dat : prior_values.txt \
	prior_position.txt \
	prior_hierarchical.txt \
	synthetic_obs_checkerboard.txt
	mkdir -p results_checkerboard_opt
	../regressionvoronois2 -i synthetic_obs_checkerboard.txt -o results_checkerboard_opt/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-T 1000 \
	-t 1000 \
	-v 10 \
	-C 500 \
	-z 10000 \
	-Z 5.0

results_checkerboard_opt/mean.txt : results_checkerboard_opt/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 0 \
	-z 0 -Z 10 \
	-T results_checkerboard_opt/stddev.txt \
	-W 90 -H 45

results_checkerboard_pt/ch.dat-000 : prior_values.txt \
	prior_position.txt \
	prior_hierarchical.txt \
	synthetic_obs_checkerboard.txt
	mkdir -p results_checkerboard_pt
	mpirun -np 8 ../regressionvoronois2pt -i synthetic_obs_checkerboard.txt -o results_checkerboard_pt/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-H prior_hierarchical.txt \
	-T 1000 \
	-t 10000 \
	-v 10 \
	-C 50 \
	-z 100 \
	-Z 5.0 \
	-c 4

results_checkerboard_pt_cont/ch.dat-000 : results_checkerboard_pt/ch.dat-000
	mkdir -p results_checkerboard_pt_cont
	mpirun -np 8 ../regressionvoronois2pt -i synthetic_obs_checkerboard.txt -o results_checkerboard_pt_cont/ \
	-I results_checkerboard_pt/ \
	-P prior_values.txt \
	-M prior_position.txt \
	-H prior_hierarchical.txt \
	-T 1000 \
	-t 10000 \
	-v 10 \
	-c 4
