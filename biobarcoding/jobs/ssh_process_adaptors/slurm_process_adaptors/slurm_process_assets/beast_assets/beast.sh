#!/bin/bash
#Acordarse de meter en ~/.beast/2.6/BEAST/templates el myTreePriors.xml y en $BEAST_DEPENDENCIES_PATH/beast/templates el myTemplate.xml
if [[ $SLURM_CLUSTER_NAME == teidehpc ]]; then
  export BEAST_DEPENDENCIES_PATH=/home/ngd_dev/beast_dependencies
elif [[ $SLURM_CLUSTER_NAME == cluster ]]; then
  export BEAST_DEPENDENCIES_PATH=/home/ngd_dev/beast_dependencies
fi


check_error()
{
  if [ $1 -ne 0 ]; then
    echo "Failure. Exit Status: $1"
    exit $1
  fi
}

export my_dir=$(pwd)
check_error $?

python3 complete_beasy_template.py "$alignments" BirthDeathModel $threads "$monophyly"
check_error $?
$BEAST_DEPENDENCIES_PATH/beast/bin/applauncher BeasyInterpreter -in beast_BirthDeathModel.bea -out beast_BirthDeathModel_tmp.xml
check_error $?
python3 estimateClockRates.py "$alignments" BirthDeathModel
check_error $?
python3 create_path_sampling.py beast_BirthDeathModel.xml
check_error $?

python3 complete_beasy_template.py "$alignments" YuleModel $threads "$monophyly"
check_error $?
$BEAST_DEPENDENCIES_PATH/beast/bin/applauncher BeasyInterpreter -in beast_YuleModel.bea -out beast_YuleModel_tmp.xml
check_error $?
python3 estimateClockRates.py "$alignments" YuleModel
check_error $?
python3 create_path_sampling.py beast_YuleModel.xml
check_error $?

mkdir coupled_mcmc coupled_mcmc/BirthDeathModel coupled_mcmc/YuleModel
check_error $?
mkdir path_sampling path_sampling/BirthDeathModel path_sampling/YuleModel
check_error $?

cd $my_dir/coupled_mcmc/BirthDeathModel
check_error $?
mv $my_dir/beast_BirthDeathModel.xml .
check_error $?
cp $my_dir/beast_loop.sh .
check_error $?
srun ./beast_loop.sh beast_BirthDeathModel.xml $threads
check_error $?
$BEAST_DEPENDENCIES_PATH/beast/bin/treeannotator -b 10 treepartition.trees BirthDeathModel.nexus
check_error $?
python3 $my_dir/nexus2newick_beast_annotations.py BirthDeathModel.nexus birthDeathModel.nexus
check_error $?

cd $my_dir/coupled_mcmc/YuleModel
check_error $?
mv $my_dir/beast_YuleModel.xml .
check_error $?
cp $my_dir/beast_loop.sh .
check_error $?
srun ./beast_loop.sh beast_YuleModel.xml $threads
check_error $?
$BEAST_DEPENDENCIES_PATH/beast/bin/treeannotator -b 10 treepartition.trees YuleModel.nexus
check_error $?
python3 $my_dir/nexus2newick_beast_annotations.py YuleModel.nexus yuleModel.nexus
check_error $?

cd $my_dir/path_sampling/BirthDeathModel
check_error $?
mv $my_dir/beast_BirthDeathModel_path_sampling.xml .
check_error $?
cp $my_dir/beast_loop.sh .
check_error $?
srun nohup bash -c "$BEAST_DEPENDENCIES_PATH/beast/bin/beast -beagle_CPU -threads $threads beast_BirthDeathModel_path_sampling.xml" > BirthDeathModel_path_sampling.out 2>&1
check_error $?

cd $my_dir/path_sampling/YuleModel
check_error $?
mv $my_dir/beast_YuleModel_path_sampling.xml .
check_error $?
cp $my_dir/beast_loop.sh .
check_error $?
srun nohup bash -c "$BEAST_DEPENDENCIES_PATH/beast/bin/beast -beagle_CPU -threads $threads beast_YuleModel_path_sampling.xml" > YuleModel_path_sampling.out 2>&1
check_error $?

cd $my_dir
check_error $?
python3 create_bayes_factor.py $my_dir/path_sampling/BirthDeathModel/BirthDeathModel_path_sampling.out $my_dir/path_sampling/YuleModel/YuleModel_path_sampling.out

