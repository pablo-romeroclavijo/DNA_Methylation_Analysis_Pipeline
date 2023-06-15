This analysis pipeline is design to process DNA methylation results originated from Illumina Next Generation Sequencing. Approximately 20.000 sequences at the time.

The programme is divided in 3 sections:

	1) Decoder: demultiplexes the data into subsamples based on barcoding of the DNA sequencing.  Data is stores in a new directory 1_Decoded_data. Runtime: 2min approx.

	2) Analyser: for every subsample, it will filter and compare the individual sequence with a provided reference. It will calculate several biologically relevant parameter and correct for technique bias. Data is stored on individual data files 2_Methylation results.

	3) Plotting: for this example very basic plots have been included. Including QC plots. Data stored in new Directory 3_Plot.


Required directories in the same directory as the script:

         + 0_Raw_data:
		1)Fastq files, can be modified in the SeqIO parameters in case of using FASTA
         + 0_Reference
		 2) Ref_uncon.txt Ñ Fasta sequence of the unconverted amplicon
		 3) Barcodes.txt Ñ file required for Decoder


Partial Scripts to generate additional sample data based on user-determined parameters can be found in SubScripts. As well as, segments of the main script.

