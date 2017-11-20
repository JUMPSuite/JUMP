#!/bin/bash
#
this_scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
dtn=`date +%Y%m%d.%N`
pwd="$(pwd)"

cd $this_scripts_dir
cd ..
SCRIPT_PARENT_DIR="$(pwd)"
cd $pwd
#
WORK_ROOT_DIR=${pwd}
#WORK_ROOT_DIR=/scratch_space/${USER}/jump
echo "SCRIPT_PARENT_DIR=${SCRIPT_PARENT_DIR}"
#
work=${WORK_ROOT_DIR}/${dtn}
touch 0
rm 0
ln -s ${dtn} 0

mkdir -p ${work}/output

cp -R $SCRIPT_PARENT_DIR/src/Class ${work}
cp -R $SCRIPT_PARENT_DIR/src/Set ${work}
cp -R $SCRIPT_PARENT_DIR/src/Spiders ${work}
cp $SCRIPT_PARENT_DIR/src/*.pl  ${work}
cp $SCRIPT_PARENT_DIR/example/jump_sj.params  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.mdx  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.prdx  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.sdx ${work}
cp ~/data/proteomics/HH_tmt10_human_jump.mzXML ${work}

#
#
#
cd ${work}
pwd
ls -al

#bsub \
#  -P proteomics-jump-${USER} \
#  -q pcgp_dev \
#  -R "rusage[mem=8000]" \
#  -eo ${work}/output/oe.log \
#  -oo ${work}/output/oo.log \
#  -cwd ${work} \
#    	${work}/jump_sj-alt.pl \
#		-p ${work}/jump_sj.params \
#		${work}/HH_tmt10_human_jump.mzXML
./jump_sj.pl \
	-p jump_sj.params \
	HH_tmt10_human_jump.mzXML		
#
#
#
