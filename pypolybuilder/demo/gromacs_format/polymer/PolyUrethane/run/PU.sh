#!/bin/bash

# Here we define some important variables.

HERE=`pwd`
PROGRAM=/home/kyam/Programs/gromacs-2018.8/install/bin/gmx_mpi
WORKDIR=${HERE}/tmp_$$

CASE=PU3x

TOPO=${WORKDIR}/${CASE}.top
MDP=${HERE}/mdp

mkdir ${WORKDIR}
cp ${HERE}/${CASE}.top ${WORKDIR}
cp ${HERE}/${CASE}.out.itp ${WORKDIR}/${CASE}.itp
cp ${HERE}/${CASE}.out.gro ${WORKDIR}/${CASE}.gro
cd ${WORKDIR}

#################### BOX ######################

${PROGRAM} editconf \
	-f ${WORKDIR}/${CASE}.gro \
	-c \
	-d 1.0 \
	-bt cubic \
	-o box.gro
	
################## SOLVATE ####################

${PROGRAM} solvate \
	-cp ${WORKDIR}/box.gro \
	-cs spc216.gro \
	-p ${TOPO} \
	-o solv.gro
	
#################### ION ######################

${PROGRAM} grompp \
	-f ${MDP}/ion.mdp \
	-c ${WORKDIR}/solv.gro \
	-p ${TOPO} \
    -maxwarn 1 \
	-o ion.tpr


echo SOL | ${PROGRAM} genion \
	-s ${WORKDIR}/ion.tpr \
	-p ${TOPO} \
	-pname NA+ \
	-nname CL- \
	-neutral \
	-o ion.gro

############ ENERGY MINIMIZATION ##############

${PROGRAM} grompp \
    -f ${MDP}/em.mdp \
    -c ${WORKDIR}/ion.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em.tpr

${PROGRAM} mdrun \
    -s em.tpr \
    -deffnm em

######### EQUILIBRATION 1: NVT  ##############

${PROGRAM} grompp \
    -f ${MDP}/nvt.mdp \
    -c ${WORKDIR}/em.gro \
    -p ${TOPO} \
    -o nvt.tpr \
    -maxwarn 2 \
    -v

${PROGRAM} mdrun \
    -s nvt.tpr \
    -deffnm nvt\
    -v

######### EQUILIBRATION 2: NPT  ##############

${PROGRAM} grompp \
    -f ${MDP}/npt.mdp \
    -c ${WORKDIR}/nvt.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o npt.tpr \
    -v

${PROGRAM} mdrun \
    -s npt.tpr \
    -deffnm npt \
    -v

########### MOLECULAR DYNAMICS ###############

${PROGRAM} grompp \
    -f ${MDP}/md.mdp \
    -c ${WORKDIR}/npt.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o md.tpr \
    -v

${PROGRAM} mdrun \
    -s md.tpr \
    -deffnm md \
    -v

#mkdir -p ${HERE}/em ${HERE}/nvt ${HERE}/npt ${HERE}/md ${HERE}/out
#mv em.* ${HERE}/em
#mv nvt.* ${HERE}/nvt
#mv npt.* ${HERE}/npt
#mv md.* ${HERE}/md
#mv *.* ${HERE}/out
#
#rm -rf ${WORKDIR}
