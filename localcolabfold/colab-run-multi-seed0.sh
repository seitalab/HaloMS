#!/bin/bash



# Please enter the your data path in XXXX, YYYY, ZZZZ.

SEED=0
CORE_FOLDER='XXX' #/export/hdd2/localcolab/human_orf/
INPUT_FOLDER=$CORE_FOLDER'input/YYY'
PROJECT_FOLDER_NAME='ZZZZ'
OUT_DIR=$CORE_FOLDER'output/'$PROJECT_FOLDER_NAME
FILE_LIST=$INPUT_FOLDER'/*'

echo $FILE_LIST;

for f in $FILE_LIST; do
    INPUT_FILE_PATH=$f
    INPUT_FILE_NAME_LIST=(${INPUT_FILE_PATH//\// })
    INPUT_FILE_NAME=${INPUT_FILE_NAME_LIST[-1]/.fasta/}
    OUT_DIR_PATH=$OUT_DIR'/'$INPUT_FILE_NAME'/'

    colabfold_batch --amber --use-gpu-relax --templates --num-recycle 3 $INPUT_FILE_PATH $OUT_DIR_PATH --random-seed 0 --num-models 1 --model-order 1 --model-type AlphaFold2-multimer-v2;      
done;


