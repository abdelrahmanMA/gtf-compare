# gtf-compare



 -   **Introduction**:

A crucial step after any task, is to evaluate the output of this task and gain insight on how efficient and accurate the process went. The task in hand is to assess the GTF output of an RNA-seq assembler by comparing it to a reference annotation and analyzing how exons, isoforms and genes are matched.

The current -and only- program that performs this task is cuffcompare from cufflinks. Cuffcompare is concerned about calculating the sensitivity and specificity on various levels including nucleotides, exons, introns, intronic chains, transcripts and locus. It can perform the analysis over several outputted GTFs from different experiments and capture the shared transcripts between them, with or without a reference annotation.

In case of provided reference annotation, cuffcompare outputs two files, .refmap and .tmap. In these two files, the program is reporting the overlap between the reference and assembly and vice versa, while stating the type of overlap. In .refmap file, the program reports each reference transcript once, and all overlapping assembled transcripts with the type of overlap as partial (c) or complete (=) match. In .tmap, the same is done but for the assembled transcripts where each transcript is reported only once with the most closely related reference transcript. The specific overlap type is then reported in this type as one of many types including complete or partial match, merely intronic overlap, generic exonic overlap, potential novel isoform with at least one splice junction shared with the reference, and many other types reported in the program documentation (1).

Another two informative file are outputted from cuffcompare which is .combined and .tracking. In .combined file, all transcripts appeared at least one time in one of the GTFs files is reported, and reported isoforms get reported only once. The .tracking file matches the transcripts between samples by satisfying specific condition such as coordinates and order of all exons, with a slack for the first and last exon.

In terms of transcript analysis, the program performs in a very informative way and provides detailed information on transcripts overlaps and overlap types. However, when it comes to further analysis or breaking down these overlaps to recognize which regions were misplaced, what boundaries where dislocated, or the overall performance of splice junctions locating, cuffcompare fails to provide this type of data. The need for further analysis performed on the exonic, junctions and gene levels to be able to give a more descriptive judge on the assemblers output and to study further reasons that could arise by examining these levels, led us to propose another program to analyze GTFs in more descriptive ways.

 - **Proposal** 
 
 **Task**: reproduce cuffCompare results with more efficiency, clarity and tabulation.

 **Description**: In order to be able to judge the efficiency of an assembler, four levels of analysis are needed:
1.  Are the exons found?
    
2.  How accurate their boundaries?
    
3.  Are the isoforms found ?
    
4.  Are the genes found?
    

In order to be able to answer each of the above questons, further analysis are required for each level where we can summarize them as follows:
 
**How to judge exons?**
1.  Ref Not found. This might happen in any of these cases:
	1.  No overlap between the ref exon and the new asm.
	2.  Minimal exonic overlap (but not “enough”)
	3.  Enough overlap but not the best match
2.  The ref exon was found with:
	1.  Imprecise boundaries (partial match)
	2.  Precise boundaries (perfect match)
3.  Novel Asm. exon. This might happen if:
	1.  No overlap with any ref exon 
	2.  Minimal exonic overlap (but not “enough”)
	
**What are the types of exons?**
1.  The first exon in the isoform with the longest 5’ end in a gene
2.  The last exon in the isoform with the longest 3’ end in a gene    
3.  Middle exon in all isoforms    
4.  Middle exon in some isoforms but the first exon in other isoforms    
5.  Middle exon junction in some isoforms but the last exon in other isoforms
 
**Note 1**: The algorithm of reference exon typing: For each reference gene, identify all its exons then collect their related metadata and define their types. For each new assembly gene, identify all its exons then collect their related metadata and define their matching pattern to the reference.

**Note 2**: The algorithm of exon matching: For each new assembled exon, find all overlapping reference exons. If it overlaps with many reference exons, assign this exon to the reference exon with longest overlap (as an absolute no of shared base pairs). If it shares the same length of sequence with multiple reference exons, choose the shortest reference exon. Now you can construct this table:

![Exons table format](https://lh3.googleusercontent.com/7hEwx5o2xxkVeZygmNCwEeQ1rBLsdmMIqiCk4weF1eyC1_E6E2bZd09QykCtiIbIDx8XKhJKZEI "Exons table format")

* Note: The same reference exon can be assigned to many Asm. exons (many of them will have partial matches)

* We can add extra step that test all the “Not Found” reference exons to define the exact subtype. Suggestion: We can make 2 versions of this table; one shows the “Novel” and “Not found” subtypes while the 2nd is concise and shows them as main category. This table will not have info about insignificant or not best matches so the last line in the above table will be removed.

  

**How to judge exonic boundaries?**
1.  Not found: if the whole reference exon was not found  
2.  Ref junction Found
	1.  Imprecise: if the ref exon is found (imprecise) but the junction is imprecise    
	2.  Precise: if the ref exon is found (precise or imprecise) and the junction is precise.    
3.  Novel junction: if the asm exon is novel  
  
**What are the types of exonic boundaries?**
1.  The start of the first exon in the isoform with the longest 5’ end in a gene    
6.  The end of the last exon in the isoform with the longest 3’ end in a gene    
7.  Exon-intron junction in all isoforms   
8.  Exon-intron junction in some isoforms but the start of the first exon in other isoforms   
9.  Exon-intron junction in some isoforms but the end of the last exon in other isoforms

**Note:** The algorithm of junction assessment: For each row in the exon output table, define the possible junctions, type of the reference junctions and do the match. This algorithm should not need to check the annotation files anymore because all calculations will depend on the exon output table

![Boundaries Table format](https://lh3.googleusercontent.com/T0e-DU3DHkFr16Kpg-X7lH35Bu2saKpfob9N-RJM5g_6m6nYZGeA_AvjYHFu_SpSWPgbOK0xotU)

**How to judge isoforms?**

 1.  Ref Not found. This might happen in any of these cases:
		1.  No exonic sequence overlap & no locus overlap    
		2.  No exonic overlap but locus overlap (e.g. intronic containment    
		3.  Minimal exonic overlap (but not “enough”)    
		4.  Enough overlap but not the best match    
 2. The ref isoform was found with:    
		1.  Imprecise intronic chain    
		2.  Precise intronic chain but wrong isoform boundaries    
		3.  Precise intronic chain and precise isoform boundaries.
    
	**Note**: Each category of the “found isoforms” can be further classified into
	 - Same strand as reference   
	 -  Opposite strand compared to the reference    
	 -  unknown strand (unknown ref, unknown asm or unknown both)    
3.  Novel isoform
	1.  Sharing insignificant sequence with reference isoform or sharing significant sequence but it is not the best one-to-one relationship    
	2.  Does not share any sequence with any reference isoform

**Note 1**: The algorism of isoform matching: For each new assembled isoform, find all overlapping reference isoforms (overlapping means share exonic sequence). If it overlaps with many reference isoforms, assign this isoform to the reference isoform with longest overlap (as an absolute no of shared base pairs). If it shares the same length of sequence with multiple reference isoforms, choose the shortest reference isoform

![Isoforms table format](https://lh3.googleusercontent.com/QGjqw5-D3S5MmgcJ-hvNzjB68KTFGTtH-yL7TbRm9mi8PDW-r7JGBCcI6RaMrCYjTbkC_XYMlBs)

* Note: The same reference isoform can be assigned to many Asm. isoforms (many of them will have partial matches)

* We can add extra step that test all the “Not Found” reference isoforms to define the exact subtype. Suggestion: We can make 2 versions of this table; one shows the “Novel” and “Not found” subtypes while the 2nd is concise and shows them as main category. This table will not have info about insignificant or not best matches so the last line in the above table will be removed.

  

**How to judge gene?**

We can use the exon and transcript output tables to generate statistics about the genes

  

**Problems**:
1.  One of the problem that could face us while implementing the code is the fact of shared isoforms between different genes. If the some exon in isoform in gene x was shared with another isoform in gene y, an in the assembler, we found the isoform of gene x. Will we report that gene y isoform exists as well?  
   **Solution**: The solution to this problem will be converting the problem from many-to-many relationship to one-to-one relationship. Instead of reporting all the isoforms overlapped between both reference and sequence, only isoform with the maximum overlapping is reported.    
2.  The second problem is the case where the best match for isoform in the reference, is not the same best match for the isoform in the assembly.  
    **Solution**: A primary solution will be to start the comparison from one file and take the best match from this file’s side.
    

**Implementation**: To implement this code, we will need to convert the information in the GTF file to bunch of features for each judgment (exons, isoforms, and genes), and for each file as well (reference and assembly). By tracing the features from the reference in the assembly, we will get a hint on the sensitivity of the assembler, while tracing the features of the assembly will give us a hint on the assembler’s specificity.

- **Flowchart**:

![Code flowchart](https://lh3.googleusercontent.com/r3wWtuleybBQNpEiDXkVa2ECZZH24p6gm2RYUdywlx52kVf1pw25SpOfZJv-_BuUb_A7KkdsiKo)
_______________________________

1.  cole-trapnell-lab.github.io/cufflinks/cuffcompare/
