#!/bin/bash

bin=../build/main
dataset=../../dataset/am

rm -f $dataset.edgepart*.8 $dataset.vertexpart*.8


$bin -p 8 -method bpart -filename $dataset -write true
$bin -p 8 -method dbh -filename $dataset -write true
$bin -p 8 -method ebv -filename $dataset -write true
$bin -p 8 -method fennel -filename $dataset -write true
$bin -p 8 -method fsm_ne -k 3 -filename $dataset -write true
$bin -p 8 -method fsm_hep -k 3 -hdf 100 -filename $dataset -write true
$bin -p 8 -method hdrf -filename $dataset -write true
$bin -p 8 -method hep -hdf 100 -filename $dataset -write true
$bin -p 8 -method hep -hdf 10 -filename $dataset -write true
$bin -p 8 -method hep -hdf 1 -filename $dataset -write true
$bin -p 8 -method hybridbl -filename $dataset -write true
$bin -p 8 -method ne -filename $dataset -write true
$bin -p 8 -method v2e_bpart -filename $dataset -write true
$bin -p 8 -method v2e_fennel -filename $dataset -write true

wc -l $dataset".edgepart.bpart.8"
wc -l $dataset".edgepart.dbh.8"
wc -l $dataset".edgepart.ebv.8"
wc -l $dataset".edgepart.fennel.8"
wc -l $dataset".edgepart.fsm_ne_k_3.8"
wc -l $dataset".edgepart.fsm_hep_k_3.8"
wc -l $dataset".edgepart.hdrf.8"
wc -l $dataset".edgepart.hep_hdf_100.8"
wc -l $dataset".edgepart.hep_hdf_10.8"
wc -l $dataset".edgepart.hep_hdf_1.8"
wc -l $dataset".edgepart.hybridbl.8"
wc -l $dataset".edgepart.ne.8"
