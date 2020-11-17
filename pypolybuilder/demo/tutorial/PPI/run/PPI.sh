#!/bin/bash

# Here we define some important variables.

HERE=/home/mayk/Documents/Labmmol/pyPolyBuilder/Tutorial/PPI/run
PROGRAM=/usr/local/gromacs-514/bin/gmx_514
WORKDIR=/tmp/mayk_$$

CASE=PPI
CARGA=0

TOPO=${WORKDIR}/${CASE}.top
MDP=${HERE}/mdp

mkdir ${WORKDIR}
cp ${HERE}/${CASE}.top ${WORKDIR}
cp ${HERE}/${CASE}.itp ${WORKDIR}
cd ${WORKDIR}

#################### BOX ######################

${PROGRAM} editconf \
	-f ${HERE}/${CASE}.gro \
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
	-o ion.tpr


${PROGRAM} genion \
	-s ${WORKDIR}/ion.tpr \
	-p ${TOPO} \
	-pname NA+ \
	-nname CL- \
	-nn ${CARGA} \
	-o ion.gro

############ ENERGY MINIMIZATION ##############

${PROGRAM} grompp \
    -f ${MDP}/em.mdp \
    -c ${WORKDIR}/ion.gro \
    -p ${TOPO} \
    -maxwarn 5 \
    -o em.tpr \
    -v

${PROGRAM} mdrun \
    -s em.tpr \
    -deffnm em \
    -v

######### EQUILIBRATION 1: NVT  ##############

${PROGRAM} grompp \
    -f ${MDP}/nvt.mdp \
    -c ${WORKDIR}/em.gro \
    -p ${TOPO} \
    -o nvt.tpr \
    -v

${PROGRAM} mdrun \
    -s nvt.tpr \
    -deffnm nvt \
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

mkdir -p ${HERE}/em ${HERE}/nvt ${HERE}/npt ${HERE}/md ${HERE}/out
mv em.* ${HERE}/em
mv nvt.* ${HERE}/nvt
mv npt.* ${HERE}/npt
mv md.* ${HERE}/md
mv *.* ${HERE}/out

rm -rf ${WORKDIR}
