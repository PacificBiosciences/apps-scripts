#! /usr/bin/bash

SUBREADSET=$1
ADAPTERS=$2
BARCODES=$3
REFERENCE=$4
OUTDIR=$5

NPROC=24
CCSPROC=8
MINBARCODESCORE=26

#Check for requires executables
RECALLADAPTERS=$(command -v recalladapters)
if [ ! -x "${RECALLADAPTERS}" ]; then
 echo 'Not found: *recalladapters* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
 exit 1
fi
LIMA=$(command -v lima)
if [ ! -x "${LIMA}" ]; then
 echo 'Not found: *lima* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
 exit 1
fi  
CCS=$(command -v ccs)
if [ ! -x "${CCS}" ]; then
 echo 'Not found: *ccs* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
 exit 1
fi
PBMM2=$(command -v pbmm2)
if [ ! -x "${PBMM2}" ]; then
 echo 'Not found: *pbmm2* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
 exit 1
fi
PARALLEL=$(command -v parallel)
if [ ! -x "${PARALLEL}" ]; then
 echo 'Not found: *gnu parallel* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
 exit 1
fi

#get movie name
MOVIENAME=$(basename ${SUBREADSET} | cut -d. -f1)

#Recall adapters
refarmDir="$(readlink -f ${OUTDIR})/refarm"
refarmPrefix="${refarmDir}/$(basename ${SUBREADSET} .subreadset.xml).refarm"
refarmSh="${refarmDir}/runRefarm.sh"

mkdir -p "${refarmDir}"
echo "${RECALLADAPTERS} \
-o ${refarmPrefix} \
-j ${NPROC} \
--adapters=$(readlink -f ${ADAPTERS}) \
--adpqc \
$(readlink -f ${SUBREADSET})" > "${refarmSh}"
bash "${refarmSh}"

#Barcoding, assume lima in path
barcodeDir="$(readlink -f ${OUTDIR})/barcode"
barcodeOut="${barcodeDir}/$(basename ${SUBREADSET} .subreadset.xml).refarm.barcoded.subreadset.xml"
barcodeSh="${barcodeDir}/runLima.sh"

mkdir -p "${barcodeDir}"
echo "${LIMA} \
--min-score ${MINBARCODESCORE} \
--split-bam-named \
--single-side \
--num-threads ${NPROC} \
${refarmPrefix}.subreadset.xml \
$(readlink -f ${BARCODES}) \
${barcodeOut}" > "${barcodeSh}"
bash "${barcodeSh}"

#CCS generation, assume gnu parallel in path
ccsDir="$(readlink -f ${OUTDIR})/ccs"
ccsSh="${ccsDir}/runCCS.sh"
ccsOut="${ccsDir}"/'$(basename {} .subreadset.xml).ccs.consensusreadset.xml'

mkdir -p "${ccsDir}"

nParallel=$(echo "$NPROC / $CCSPROC" | bc)
if [ "${nParallel}" -lt 1 ]; then
 nParallel=1
 nproc="${NPROC}"
else
 nproc=${CCSPROC}
fi

echo "${PARALLEL} \
-j ${nParallel} \
${CCS} \
--numThreads ${nproc} \
--minPasses 3 \
--minPredictedAccuracy 0.99 \
{} '${ccsOut}' \
::: ${barcodeDir}/*.subreadset.xml" > "${ccsSh}"
bash "${ccsSh}"

#Align ccs
alignDir="$(readlink -f ${OUTDIR})/align"
alignSh="${alignDir}/runAlign.sh"
bcList=$(ls ${ccsDir}/*xml | awk -F. '{print $(NF-3)}' | paste -s)

mkdir -p "${alignDir}"

echo "${PARALLEL} \
-j ${nParallel} \
${PBMM2} align \
--sort \
--alignment-threads ${nproc} \
-r 2000 \
-A 2 \
-B 5 \
-z 400 \
-Z 50 \
-o 5 \
-O 56 \
-e 4 \
-E 0 \
-L 0.2 \
-c 0 \
--sample {} \
$(readlink -f $REFERENCE) \
${ccsDir}/*{}.ccs.consensusreadset.xml \
${alignDir}/${MOVIENAME}.refarm.barcoded.aligned.{}.consensusalignmentset.xml \
::: ${bcList}" > "${alignSh}"
bash "${alignSh}"

