SUBREADSET=$1
ADAPTERS=$2
BARCODES=$3
REFERENCE=$4
OUTDIR=$5

#Recall tool, not currently available for customers
RECALLADAPTERS=/mnt/software/p/ppa/develop-61955/bam2bam
NPROC=48
CCSPROC=8
MINBARCODESCORE=26

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
echo "$(which lima) \
--same \
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
ccsOut="${ccsDir}"/'$(basename {} .subreadset.xml).consensusreadset.xml'

mkdir -p "${ccsDir}"

nParallel=$(echo "$NPROC / $CCSPROC" | bc)
if [ "${nParallel}" -lt 1 ]; then
 nParallel=1
 nproc="${NPROC}"
else
 nproc=${CCSPROC}
fi

echo "$(which parallel) \
-j ${nParallel} \
$(which ccs) \
--numThreads ${nproc} \
{} '${ccsOut}' \
::: ${barcodeDir}/*.subreadset.xml" > "${ccsSh}"
bash "${ccsSh}"

#Align ccs
alignDir="$(readlink -f ${OUTDIR})/align"
alignSh="${alignDir}/runAlign.sh"
bcList=$(ls ${ccsDir}/*xml | awk -F. '{print $(NF-2)}' | paste -s)

mkdir -p "${alignDir}"

echo "$(which parallel) \
-j ${nParallel} \
$(which pbmm2) align \
--sort \
--alignment-threads ${nproc} \
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
--sample {} \
$(readlink -f $REFERENCE) \
${ccsDir}/*{}.consensusreadset.xml \
${alignDir}/${MOVIENAME}.refarm.barcoded.{}.consensusalignmentset.xml \
::: ${bcList}" > "${alignSh}"
bash "${alignSh}"

