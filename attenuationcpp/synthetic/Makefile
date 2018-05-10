
#
# Q is between 0 and 1000 roughly speaking hence log(Q) ranges from -3 to 7
# with the lower bound being fairly arbitrary (0 would probably be reasonable
#


synthetic_constant_0.1.txt : ../data/coreattenuation_converted.txt
	../mksynthetictomography -i $< -o $@ \
	-m Constant \
	-s 7.0 \
	-n 0.1 \
	-I synthetic_constant_0.1.img

synthetic_checkerboard_0.1.txt : ../data/coreattenuation_converted.txt
	../mksynthetictomography -i $< -o $@ \
	-m CheckerBoard \
	-s 7.0 \
	-n 0.1 \
	-I synthetic_checkerboard_0.1.img

value_prior.txt :
	echo "prior Uniform" > $@
	echo "0.0 8.0" >> $@
	echo "proposal Gaussian" >> $@
	echo "0.05" >> $@

position_prior.txt :
	echo "sphericalprior UniformSpherical" > $@
	echo "sphericalproposal VonMisesSpherical" >> $@
	echo "2.0e2" >> $@

results_constant_opt/ch.dat : value_prior.txt \
	position_prior.txt \
	synthetic_constant_0.1.txt
	mkdir -p results_constant_opt
	../attenuationvoronois2 -i synthetic_constant_0.1.txt -o results_constant_opt/ \
	-P value_prior.txt \
	-M position_prior.txt \
	-T 100 \
	-t 1000 \
	-v 10 \
	-C 10 \
	-z 10000 \
	-Z 5.0

results_constant_opt/mean.txt : results_constant_opt/ch.dat
	../../postS2Voronoi_mean -i $< -o $@ \
	-t 10 -s 0 \
	-z 0 -Z 10 \
	-T results_constant_opt/stddev.txt \
	-W 90 -H 45

results_checkerboard_opt/ch.dat : value_prior.txt \
	position_prior.txt \
	synthetic_checkerboard_0.1.txt
	mkdir -p results_checkerboard_opt
	../attenuationvoronois2 -i synthetic_checkerboard_0.1.txt -o results_checkerboard_opt/ \
	-P value_prior.txt \
	-M position_prior.txt \
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




