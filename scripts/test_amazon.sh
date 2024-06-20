#!/bin/bash

bin=../build/main
dataset=../../dataset/am

#############################
#### TEST ONEFILE OUTPUT ####
#############################

$bin -p 8 -method bpart -filename $dataset -write onefile
$bin -p 8 -method dbh -filename $dataset -write onefile
$bin -p 8 -method ebv -filename $dataset -write onefile
$bin -p 8 -method fennel -filename $dataset -write onefile
$bin -p 8 -method fsm_ne -k 3 -filename $dataset -write onefile
$bin -p 8 -method fsm_hep -k 3 -hdf 100 -filename $dataset -write onefile
$bin -p 8 -method hdrf -filename $dataset -write onefile
$bin -p 8 -method hep -hdf 100 -filename $dataset -write onefile
$bin -p 8 -method hep -hdf 10 -filename $dataset -write onefile
$bin -p 8 -method hep -hdf 1 -filename $dataset -write onefile
$bin -p 8 -method hybridbl -filename $dataset -write onefile
$bin -p 8 -method ne -filename $dataset -write onefile
$bin -p 8 -method v2e_bpart -filename $dataset -write onefile
$bin -p 8 -method v2e_fennel -filename $dataset -write onefile

wc -l $dataset.edgepart.bpart.8
wc -l $dataset.edgepart.dbh.8
wc -l $dataset.edgepart.ebv.8
wc -l $dataset.edgepart.fennel.8
wc -l $dataset.edgepart.fsm_ne_k_3.8
wc -l $dataset.edgepart.fsm_hep_k_3.8
wc -l $dataset.edgepart.hdrf.8
wc -l $dataset.edgepart.hep_hdf_100.8
wc -l $dataset.edgepart.hep_hdf_10.8
wc -l $dataset.edgepart.hep_hdf_1.8
wc -l $dataset.edgepart.hybridbl.8
wc -l $dataset.edgepart.ne.8

rm -f $dataset.edgepart*.8 $dataset.vertexpart*.8


#############################
### TEST MULTIFILE OUTPUT ###
#############################

$bin -p 8 -method dbh -filename $dataset -write multifile
$bin -p 8 -method ebv -filename $dataset -write multifile
$bin -p 8 -method fsm_ne -k 3 -filename $dataset -write multifile
$bin -p 8 -method fsm_hep -k 3 -hdf 100 -filename $dataset -write multifile
$bin -p 8 -method hdrf -filename $dataset -write multifile
$bin -p 8 -method hep -hdf 100 -filename $dataset -write multifile
$bin -p 8 -method hep -hdf 10 -filename $dataset -write multifile
$bin -p 8 -method hep -hdf 1 -filename $dataset -write multifile
$bin -p 8 -method hybridbl -filename $dataset -write multifile
$bin -p 8 -method ne -filename $dataset -write multifile

wc -l $dataset.edgepart.dbh.8.part*
wc -l $dataset.edgepart.ebv.8.part*
wc -l $dataset.edgepart.fennel.8.part*
wc -l $dataset.edgepart.fsm_ne_k_3.8.part*
wc -l $dataset.edgepart.fsm_hep_k_3.8.part*
wc -l $dataset.edgepart.hdrf.8.part*
wc -l $dataset.edgepart.hep_hdf_100.8.part*
wc -l $dataset.edgepart.hep_hdf_10.8.part*
wc -l $dataset.edgepart.hep_hdf_1.8.part*
wc -l $dataset.edgepart.hybridbl.8.part*
wc -l $dataset.edgepart.ne.8.part*

rm -f $dataset.edgepart*.8.part*