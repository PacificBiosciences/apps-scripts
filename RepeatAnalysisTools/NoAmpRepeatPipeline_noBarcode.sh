SUBREADSET=$1
ADAPTERS=$2
REFERENCE=$3
OUTDIR=$4

#Recall tool, not currently available for customers
NPROC=24

#Check for requires executables
RECALLADAPTERS=$(command -v recalladapters)
if [ ! -x "${RECALLADAPTERS}" ]; then
 echo 'Not found: *recalladapters* is needed to run this script. See https://github.com/PacificBiosciences/pbbioconda'
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

#CCS generation
ccsDir="$(readlink -f ${OUTDIR})/ccs"
ccsSh="${ccsDir}/runCCS.sh"
ccsOut="${ccsDir}/${MOVIENAME}.consensusreadset.xml"

mkdir -p "${ccsDir}"

echo "${CCS} \
--numThreads ${NPROC} \
--minPasses 3 \
--minPredictedAccuracy 0.99 \
${refarmPrefix}.subreadset.xml \
${ccsOut}" > "${ccsSh}"
bash "${ccsSh}"

#Align ccs
alignDir="$(readlink -f ${OUTDIR})/align"
alignOut="${alignDir}/${MOVIENAME}.consensusalignmentset.xml"
alignSh="${alignDir}/runAlign.sh"

mkdir -p "${alignDir}"

echo "${PBMM2} align \
--sort \
--alignment-threads ${NPROC} \
-r 2k \
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
$(readlink -f $REFERENCE) \
${ccsOut} \
${alignDir}/${MOVIENAME}.refarm.consensusalignmentset.xml" > "${alignSh}"
bash "${alignSh}"
