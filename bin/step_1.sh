#!/bin/bash
#
this_work_dir=$1
dtn=$2
#

#
this_scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#

cd $this_scripts_dir
cd ..
cd $this_work_dir

#
work=${this_work_dir}/${dtn}
touch 0
rm 0
ln -s ${dtn} 0

#mkdir -p ${work}/output

cp -R $this_work_dir/src/Class ${work}
cp -R $this_work_dir/src/Set ${work}
cp -R $this_work_dir/src/Spiders ${work}
cp $this_work_dir/src/*.pl  ${work}
cp $this_work_dir/example/jump_sj.params  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.mdx  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.prdx  ${work}
cp ~/data/proteomics/human_ft_mc2_c0_TMT_K229.fasta.sdx ${work}
cp ~/data/proteomics/HH_tmt10_human_jump.mzXML ${work}

#
#
#
cd ${work}
pwd
#ls -al

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
