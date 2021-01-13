# Semi automated pipeline

1. run_prepare.sh <br />
      -runs codons.py from 'p53' dir, change this.<br />
      -Also runs Prune_nuc_and_protein.py (basedir)<br />
2. run_treemmer.sh (run from Treemmer dir or modify) <br />
3. Parse_Treemmer_Verbose2.ipynb <br />
4. CreateSubSampledUnalignedFasta.ipynb <br />
5. run_subsampled.sh<br />
      -runs Prune_nuc_and_protein.py (basedir) <br />
