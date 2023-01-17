case="17-08"
refname="hg38"
PWD="/confidential/FamilyR13/DATA/10x/sv_compare/results/${case}_${refname}"

echo $PWD

step5=false
step6=false
step7=true

if $step5; then
	mkdir -p ensemble
	python3 /confidential/tGenVar/scripts/tGenVar/dicast/00_prepare_ensemble.py ${case} ${refname} ${PWD}/${case}.decast_input.csv ${PWD}/ensemble/${case}_${refname}.SVs.raw.tsv /confidential/FamilyR13/DATA/10x/sv_compare/params_ensemble.json
fi

if $step6; then
    python3 /confidential/tGenVar/scripts/tGenVar/dicast/01_collect_reference.py ${case} ${refname} ${PWD}/ensemble/${case}_${refname}.SVs.raw.tsv ${PWD}/ensemble/${case}_${refname}.SVs.ref.tsv /confidential/FamilyR13/DATA/10x/sv_compare/params_ensemble.json
fi

if $step7; then
    for CHR in {1..22} X Y
    do
        python3 /confidential/tGenVar/scripts/tGenVar/dicast/01_collect_alignment.py ${case} ${refname} ${PWD}/ensemble/${case}_${refname}.SVs.raw.tsv ${PWD}/ensemble/${case}_${refname}.SVs.aln.ill.chr${CHR}.tsv /confidential/FamilyR13/DATA/10x/sv_compare/params_ensemble.json chr${CHR} &
    done

	trap 'kill 0' INT   # make ^C work
	status=0            # exit status of this script, assume okay
	while true; do
		wait -n                     # wait for any child
		sts=$?                      # capture exit status of wait
		(($sts == 127)) && break    # if 127, no more children
		(($sts)) && status=1        # otherwise exit status of child. if bad, propagate
	done
	
fi