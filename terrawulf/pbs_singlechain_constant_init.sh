#!/bin/bash

##########################
#                        #
#   The PBS directives   #
#                        #
##########################

#          Set the name of the job (up to 15 characters, 
#          no blank spaces, start with alphanumeric character)

#PBS -N gvs2const

#          Specify the number of nodes requested,the type of nodes and 
#          the number of processors per node.
#          In this example, we request 1 nodes of ethernet type with 
#          4 processors on each node.  

#PBS -l nodes=1:t3:ppn=12

## Alternative: PBS -l nodes=2:t3:ppn=12

#          By default, the standard output and error streams are sent
#          to files in the current working directory with names:
#              job_name.osequence_number  <-  output stream
#              job_name.esequence_number  <-  error stream
#          where job_name is the name of the job and sequence_number 
#          is the job number assigned when the job is submitted.
#          Use the directives below to change the files to which the
#          standard output and error streams are sent.

# REPLACE $USER by you username


#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output. 

#PBS -j oe

#          Specify the maximum wall clock time. The wall
#          clock time should take possible queue waiting time into
#          account.  Format:   hhhh:mm:ss   hours:minutes:seconds
#          Be sure to specify a reasonable value here.
#          If the job does not finish by the time reached,
#          the job is terminated.

#PBS -l walltime=48:00:00



#          TORQUE can send informative email messages to you about the
#          status of your job.  Specify a string which consists of
#          either the single character "n" (no mail), or one or more
#          of the characters "a" (send mail when job is aborted),
#          "b" (send mail when job begins), and "e" (send mail when
#          job terminates).  The default is "a" if not specified.
#          You should also specify the email address to which the
#          message should be send via the -M option.

#  #PBS -m abe
#  #PBS -m ae

#  #PBS -m n
#PBS -M $USER@rses.anu.edu.au



##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

NCPU=`wc -l < $PBS_NODEFILE`
echo ------------------------------------------------------
echo ' This job is allocated on '${NCPU}' cpu(s)'
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------


############################################################
#                                                          #
#    Execute the run.  Do not run in the background.       #
#    pbsdsh $PROGRAM_NAME
#                                                          #
############################################################

export MPD_CON_EXT=$PBS_JOBID

# This needs to change to point to the base directory of the software
export GENERALBASE=/home/rhys/tools/GeneralVoronoiS2

export GENERALINVERT_PT=$GENERALBASE/generalregressioncpp/regressionvoronois2pt

# These are provided as examples only and don't exist and either need to be created or
# lines edited below to point to correct files
export GENERAL_OBS=$GENERALBASE/terrawulf/observations.txt
export PRIORFILE=$GENERALBASE/terrawulf/prior.txt
export POSITIONPRIORFILE=$GENERALBASE/terrawulf/position_prior.txt
export HIERARCHICALPRIORFILE=$GENERALBASE/terrawulf/synthetic_hierarchical_prior.txt
export RESULTS_DIR=$GENERALBASE/terrawulf/constant_mpi_0/

export GENERALARGS="-i $GENERAL_OBS -o $RESULTS_DIR -P $PRIORFILE -M $POSITIONPRIORFILE -t 1000000 -v 1000 -T 200 -l 1.0 -H $HIERARCHICALPRIORFILE"

mkdir -p $RESULTS_DIR
echo $GENERALARGS > $RESULTS_DIR/run0.args

mpirun -np $NCPU -machinefile $PBS_NODEFILE $GENERALINVERT_PT $GENERALARGS 

exit
