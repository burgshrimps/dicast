# `dicast check` drop-reason glossary

Every record `dicast check` (and, internally, `dicast call`/`multi`) reads
from an input VCF is either **kept** (counted under `DEL`/`DUP`/`INS`) or
**dropped** for exactly one of the reasons below. `notes` counters below
that are informational only -- they apply alongside a *kept* record and
never drop it.

| dropped reason | meaning | what to do |
|---|---|---|
| `contig` | the record's contig could not be mapped onto any name in `--fai` (identity, `+chr`, or `-chr` were all tried and none matched) | rename contigs with `bcftools annotate --rename-chrs`, or double-check `--fai` is the right reference |
| `no_svtype` | no `INFO/SVTYPE`, no `INFO/SIMPLE_TYPE`, no recognizable symbolic or breakend ALT, and REF/ALT gave no usable sequence to compare (missing ALT, an unrecognized symbolic ALT, or the `*` overlapping-deletion allele) | the caller's output doesn't carry any of the SV-type signals dicast looks for; check the caller's documentation for how it encodes SV type |
| `INV` | inversion call | out of scope for dicast's model (DEL/DUP/INS only); not fixable, expected to be dropped |
| `BND` | breakend call with no `INFO/SIMPLE_TYPE` override (raw gridss, delly `TRA`, or any bracket-notation ALT with plain `SVTYPE=BND`) | dicast does not pair breakends itself -- see "Fixing inputs" in `SKILL.md` (gridss `simple-event-annotation.R`; for delly, use its DEL/DUP/INS calls instead) |
| `CNV` | copy-number call with no more specific type | out of scope; not fixable |
| `other_svtype` | a recognized SV type outside `{DEL, DUP, INS}` (e.g. `DUP:TANDEM`'s stripped form is `DUP` and is fine, but something like `IDUP`/custom caller-specific types are not) | out of scope for dicast's model; not fixable |
| `small_variant` | a sequence-resolved REF/ALT pair with `\|len(ALT) - len(REF)\| < 50 bp` (dicast's SV size floor) | expected for indel-sized calls; not an SV by dicast's contract |
| `missing_end` | a DEL/DUP record with no derivable end position (`rec.stop` unusable and no `INFO/SVLEN`) | the caller omitted both `END` and `SVLEN`; check whether a less-filtered version of the caller's output restores one of them |
| `bnd_mate` | a `SIMPLE_TYPE`-annotated breakend record that is the *upper* breakend of its pair (its mate is on a different contig, or at an earlier position) -- only the lower breakend of each pair is kept, to avoid double-counting the event | expected; the paired lower-breakend record carries the same event |
| `absent_in_sample` | in a multi-sample VCF, the chosen `--sample` column's `GT` has no alt allele (`> 0`) at this record | expected for cohort VCFs where not every sample carries every call; not a data problem |

## notes (informational, not drops)

| note | meaning |
|---|---|
| `ins_len_unknown` | an INS record was kept but its length could not be derived from `INFO/SVLEN` or a sequence-resolved ALT, so `sv_len` is `NaN`. dicast's XGBoost models tolerate `NaN` features, so this is not an error, just a coverage gap for that record. |
| `no_gt` | a record was kept but the chosen sample column had no `GT` field at all (a `./.` genotype does *not* count -- that is a present-but-missing GT, e.g. raw lumpy output); genotype defaults to `(None, None)`. |
