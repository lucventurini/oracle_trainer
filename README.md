This is a collection of scripts to perform the training of TwinScan/SNAP/Augustus/GeneID.

Required:

augustus	(set AUGUSTUS_CONFIG_PATH accordingly)
N-scan and IPE	(http://mblab.wustl.edu/software.html)
EVAL   	   	(http://mblab.wustl.edu/software.html; set the parameter "EVAL_GTF" to point to the directory containing its scripts)
gtf2gff3	(http://www.sequenceontology.org/software/GAL_Code/GAL_0.2.2_stand_alone_scripts.tar.gz)


As input, prepare_dataset.py takes the proteins, the GTF, and the masked and unmasked genome. The GTF *must* contain the
start/stop codon features; they can be added with refeature from PHAST (http://compgen.cshl.edu/phast/index.php).