# gtf-compare



 -   **Introduction**:

A crucial step after any task, is to evaluate the output of this task and gain insight on how efficient and accurate the process went. The task in hand is to assess the GTF output of an RNA-seq assembler by comparing it to a reference annotation and analyzing how exons, isoforms and genes are matched.

The current -and only- program that performs this task is cuffcompare from cufflinks. Cuffcompare is concerned about calculating the sensitivity and specificity on various levels including nucleotides, exons, introns, intronic chains, transcripts and locus. It can perform the analysis over several outputted GTFs from different experiments and capture the shared transcripts between them, with or without a reference annotation.

In case of provided reference annotation, cuffcompare outputs two files, .refmap and .tmap. In these two files, the program is reporting the overlap between the reference and assembly and vice versa, while stating the type of overlap. In .refmap file, the program reports each reference transcript once, and all overlapping assembled transcripts with the type of overlap as partial (c) or complete (=) match. In .tmap, the same is done but for the assembled transcripts where each transcript is reported only once with the most closely related reference transcript. The specific overlap type is then reported in this type as one of many types including complete or partial match, merely intronic overlap, generic exonic overlap, potential novel isoform with at least one splice junction shared with the reference, and many other types reported in the program documentation (1).

Another two informative file are outputted from cuffcompare which is .combined and .tracking. In .combined file, all transcripts appeared at least one time in one of the GTFs files is reported, and reported isoforms get reported only once. The .tracking file matches the transcripts between samples by satisfying specific condition such as coordinates and order of all exons, with a slack for the first and last exon.

In terms of transcript analysis, the program performs in a very informative way and provides detailed information on transcripts overlaps and overlap types. However, when it comes to further analysis or breaking down these overlaps to recognize which regions were misplaced, what boundaries where dislocated, or the overall performance of splice junctions locating, cuffcompare fails to provide this type of data. The need for further analysis performed on the exonic, junctions and gene levels to be able to give a more descriptive judge on the assemblers output and to study further reasons that could arise by examining these levels, led us to propose another program to analyze GTFs in more descriptive ways.

 - **Proposal** 
 
 **Task**: 
1.  Extend the Cuffcompare functionality from comparing isoforms (Transfrags) to compare exons, exon-exon junctions, and genes
2.  Provide better matrices to compare the differences between annotations (This should produce better summary stats than current Cuffcompare)
3.  Gene specific or interval specific reports 

**Before we start:** 
-------------------- 
**Level of restriction**

Dealing with stand issue in GTF is always tricky. I can think of 3 scenarios: 
1.  Strict: 3 lists for positive, negative and unknown strands
2.  Intermediate: 2 lists for all possible positive and negative strands. Add the intervals of unknown strand to both lists
3.  Lenient: 1 list for all intervals 
I think we should implement the “Intermediate” level for now. However we have to consider the effect of the stand on the gene annotation. Remeber that the 5' and 3' ends of an isofomrs on the +ve strand are ooposite to those on the -ve strand  

**GTF as a tree structure:**

Each gene has a unique id and is associated with a no of isoforms. Each isoform in a gene has a unique id and is associated with a no of exons.  Each exon in an isoform has a unique id. However we typically have multiple identical exons in the different isoforms so we end up with having this same exonic stretch with multiple ids. Also theoretically we might have the same isoforms with multiple ids because it belongs to multiple genes (This is not normal but can happen in crazy assembly or merge).
Can we have exons or transcripts without genes? I saw some cases like this. How does the parser behave? Should we ignore for now?    

**Types of exons**
1.  The first exon in the isoform 
	1.  With the longest 5’ end in a gene
2.  The last exon in the isoform 
	1.  With the longest 3’ end in a gene    
3.  Middle exon in all isoforms    
4.  Middle exon in some isoforms but the first exon in other isoforms    
5.  Middle exon in some isoforms but the last exon in other isoforms

**Types of exonic boundaries?**
1.  The start of the first exon in the isoform (it is not a splice junction) 
	1.  With the longest 5’ end in a gene  ==> this will be the start of the gene  
2.  The end of the last exon in the isoform (it is not a splice junction) 
	1.  With the longest 3’ end in a gene ==> this will be the end of the gene     
3.  Exon-intron junction in all isoforms   
4.  Exon-intron junction in some isoforms but the start of the first exon in other isoforms   
5.  Exon-intron junction in some isoforms but the end of the last exon in other isoforms

 **Compenets of the problem**: To judge the efficiency of an assembler, four levels of analysis are needed:
1.  Are the genes found?
2.  Are the isoforms found? How accurate their starnd, intronic chains, and boundaries?
3.  Are the exons found? How accurate their boundaries?

In order to be able to answer each of the above questons, further analysis are required for each level where we can summarize them as follows:

**How to judge gene/locus?**

We have two level:
1. Locus overlap 
2. Stats of matching for its transcripts and exons  

**How to judge isoforms?**

1.  Precise intronic chain and precise isoform boundaries.
2.  Precise intronic chain but wrong isoform boundaries.
3.  Query is contained (The query is missing junctions on either or both ends but the remaining > 0 & precise)
4.  Exonic overlap with imprecise intronic chain (but not Query containment). ==> See the note  
	- we can make subcategories here: e.g. 
	1.  Reference is contained (The query has more junctions on either or both ends but the remaining > 0 & precise)
	2.  check Cuffcompre for other options
5.  No exonic overlap but overlap with reference on the opposite strand
6.  No exonic overlap but locus overlap ==> intronic containment, novel 5' expression, or novel 3' expression	
7.  No exonic sequence overlap & no locus overlap    
 **Note**: This type require extra-output table that has:
- no of junc in query/ no of matching junc/ no of junc in ref
- no of shared exonic bases/ no of intronic bases/ no of 5' novel bases/ no of 3' novel bases
 
**How to judge exons?**
1.  Overlap with one exon with precise boundaries (perfect match)
2.  Overlap with one exon with imprecise boundaries (partial match: containment or overlapping): count base diff
3.  More than one-to-one relationship: count base diff
	1.  Reference exon overlaps with multiple query exons 
	2.  Query exon overlaps with multiple reference exons
	3.  Many reference exons overlap with many query exons
4.  No overlap between the ref exon and the new asm
**Note** It could be useful to judge the parent isoforms: What is the best relationship between their possible parents isoforms?

**How to judge exonic boundaries?**
Each junction is defined by the last base position of the 5' exon & the first base position of the 3' exon
1.  Percise known match: The junction exactly match a reference one
2.  Percise unknown match: the junction match a start and end of known exons but the junction itself is not in the ref
3.  Imprecise: The junction is half known
4.  Novel: both ends are novel 


**Algorithm** 
1.  Parse genes in the GTFs (or the DBs) into a list sorted by co-ordinates (note that gene id is unique but gene name does not need to)
2.  Link each locus with a list of its isoforms (an array of arrays for the isoform ids)
3.  Convert genes into a list of intervals per chromosome. * See the note about the list design. If the gene is not defined in the GTF and we can not query the DB for its boundaries, define the gene interval that starts by the beginning of the most 5’ isoform and ends by the end of the most 3’ isoform
	* Is it better to have separate lists per strand? If yes the gene should take the strand sign of its isoform with the highest no of splice junctions and defined strand.
4.  Find overlapping intervals:
	1.  For a one-way analysis: Find all reference intervals overlapping with each query interval
	2.  For a two-way analysis: Find the opposite as well. Do we need to consider a super locus that include all overlapping loci?
5.  Match exons:
	1.  Make a map between all exon ids in the tested interval and their co-ordinates. 
	2.  Create a unique list of exonic co-ordinates (now this is the actual list of unique exons in the gene)
	3.  Match the exons sequenctially 
6.  Parse isoforms in the GTFs (or the DBs) into 2 lists (for the 2 strands) sorted by co-ordinates
7.  Link each isoform with a list of its exon intervals (These are the unique inetrvals)
8.  Match Junctions
	1.  Make a list of all junctions in all isoform
	2.  merge duplicates 
	3.  Match query versus reference 
9. Match isoforms:  
	1. Compare each query and reference isoforms using junction lists
	2. Compare each query and reference isoforms using unique exon inetrval lists
	* Another suggestion
	1.  Convert each isoform into bed 
	2.  Calc the bed intersection of each two
10.  Locus stats

List design: I can think of 2 designs to enable efficient search:
•  Random access design: Each list can be presented as 2 arrays S and E (all structures are sorted): One for the starts and one for the ends. Now you can check for any point (x) if it falls in any interval by finding the first item (i) in E array where Ei >= x then confirm that x <= Si
I think this will be best fit for one-way search with the query items much smaller than the reference
•  Sequential access design 


Notes:
- Edge case: If the gene is not defined in the GTF, all the overlapping isoforms on the same strand will be considered one locus (un-stranded isoforms will be included with both strands). Overlapping means sharing exonic sequence. 


**Note:**: 

-  The algorithm of exon matching: For each new assembled exon, find all overlapping reference exons. If it overlaps with many reference exons, assign this exon to the reference exon with longest overlap (as an absolute no of shared base pairs). If it shares the same length of sequence with multiple reference exons, choose the shortest reference exon. Now you can construct this table:
-  The same reference isoform or exon can be assigned to many Asm. exons (many of them will have partial matches)
-  The algorithm of junction assessment: For each row in the exon output table, define the possible junctions, type of the reference junctions and do the match. This algorithm should not need to check the annotation files anymore because all calculations will depend on the exon output table
-  In the two-way mode, the “found isoforms” can be additionally labeled by its rank match aganist the reference 
  

**Problems**:
1.  One of the problem that could face us while implementing the code is the fact of shared isoforms between different genes. If the some exon in isoform in gene x was shared with another isoform in gene y, an in the assembler, we found the isoform of gene x. Will we report that gene y isoform exists as well?  
   **Solution**: The solution to this problem will be converting the problem from many-to-many relationship to one-to-one relationship. Instead of reporting all the isoforms overlapped between both reference and sequence, only isoform with the maximum overlapping is reported.    
2.  The second problem is the case where the best match for isoform in the reference, is not the same best match for the isoform in the assembly.  
    **Solution**: A primary solution will be to start the comparison from one file and take the best match from this file’s side.
    

**Implementation**: To implement this code, we will need to convert the information in the GTF file to bunch of features for each judgment (exons, isoforms, and genes), and for each file as well (reference and assembly). By tracing the features from the reference in the assembly, we will get a hint on the sensitivity of the assembler, while tracing the features of the assembly will give us a hint on the assembler’s specificity.



1.  cole-trapnell-lab.github.io/cufflinks/cuffcompare/
