MODEL=${1}

for SVTYPE in DEL INS INV DUP
do
    python3 dicast.py train ${SVTYPE} params/tgenvar/params_train.json params/model/params_model.json ${MODEL} --cur
done