# V1_TM
## Task
Reproduce cuffCompare results with more efficiency, clarity and tabulation.

## Description

> In order to be able to judge the efficiency of an assembler, four
> levels of analysis are needed, where each level requires further
> analysis.

1.  Are the exons found?
2.  How accurate their boundaries?
3.  Are the isoforms found ?
4.  Are the genes found?
## How to judge exons?

1.  **Ref Not found. This might happen in any of these cases:**
    a.  No overlap between the ref exon and the new asm.
    b.  Minimal exonic overlap (but not “enough”)
    c.  Enough overlap but not the best match

2.  **The ref exon was found with:**

    a.  Imprecise boundaries (partial match)
    b.  Precise boundaries (perfect match)
  
3.  **Novel Asm. exon. This might happen if:**

    a.  No overlap with any ref exon
    b.  Minimal exonic overlap (but not “enough”)

## What are the types of exons?

1.  The first exon in the isoform with the longest 5’ end in a gene
2.  The last exon in the isoform with the longest 3’ end in a gene
3.  Middle exon in all isoforms
4.  Middle exon in some isoforms but the first exon in other isoforms
5.  Middle exon junction in some isoforms but the last exon in other isoforms

***Note 1: The algorithm of reference exon typing:*** For each reference gene, identify all its exons then collect their related metadata and define their types. For each new assembly gene, identify all its exons then collect their related metadata and define their matching pattern to the reference.

***Note 2: The algorithm of exon matching:*** For each new assembled exon, find all overlapping reference exons. If it overlaps with many reference exons, assign this exon to the reference exon with longest overlap (as an absolute no of shared base pairs). If it shares the same length of sequence with multiple reference exons, choose the shortest reference exon. Now you can construct this table:


<table>
  <tr>
    <th colspan="4">Reference<br>Exons</th>
    <th colspan="4">New Asm. Exons</th>
  </tr>
  <tr>
    <th>Position<br></th>
    <th>Gene<br></th>
    <th>Isoform</th>
    <th>Type</th>
    <th>Position</th>
    <th>Gene</th>
    <th>Isoform</th>
    <th>Type</th>
  </tr>
  <tr>
    <td>Ref. Exon<br></td>
    <td>G1</td>
    <td>Iso1,Iso2</td>
    <td>Type 1-5</td>
    <td>--</td>
    <td>--</td>
    <td>--</td>
    <td>Not Found</td>
  </tr>
  <tr>
    <td>Ref. Exon</td>
    <td>G1</td>
    <td>Iso2</td>
    <td>Type 1-5</td>
    <td></td>
    <td>NG1</td>
    <td>NIso1</td>
    <td>Partial or Perfect</td>
  </tr>
  <tr>
    <td>---</td>
    <td></td>
    <td></td>
    <td>--</td>
    <td></td>
    <td>NG1</td>
    <td>NIso2</td>
    <td>Novel-A</td>
  </tr>
  <tr>
    <td>Ref. Exon</td>
    <td>G2</td>
    <td>Iso1</td>
    <td>Type 1-5</td>
    <td></td>
    <td>NG1</td>
    <td>NIso2</td>
    <td>Novel-B</td>
  </tr>
</table>


\* Note: The same reference exon can be assigned to many Asm. exons (many of them will have partial matches)

\*\* We can add extra step that test all the “Not Found” reference exons to define the exact subtype. Suggestion: We can make 2 versions of this table; one shows the “Novel” and “Not found” subtypes while the 2nd is concise and shows them as main category. This table will not have info about insignificant or not best matches so the last line in the above table will be removed.

## How to judge exonic boundaries?

1.  **Not found**: if the whole reference exon was not found
2.  **Ref junction Found**
    a.  Imprecise: if the ref exon is found (imprecise) but the junction is imprecise
    b.  Precise: if the ref exon is found (precise or imprecise) and the junction is precise.

3.  **Novel junction:** if the asm exon is novel

## What are the types of exonic boundaries?

1.  The start of the first exon in the isoform with the longest 5’ end in a gene
2.  The end of the last exon in the isoform with the longest 3’ end in a gene
3.  Exon-intron junction in all isoforms
4.  Exon-intron junction in some isoforms but the start of the first exon in other isoforms
5.  Exon-intron junction in some isoforms but the end of the last exon in other isoforms

> Note: The algorithm of junction assessment: For each row in the exon
> output table, define the possible junctions, type of the reference
> junctions and do the match. This algorithm should not need to check
> the annotation files anymore because all calculations will depend on
> the exon output table


<table>
  <tr>
    <th colspan="2">Reference Assembly</th>
    <th colspan="2">New Assembly</th>
  </tr>
  <tr>
    <th>Junction Position</th>
    <th>Junction Type</th>
    <th>Junction Position</th>
    <th>Relation</th>
  </tr>
  <tr>
    <td>Reference Junction</td>
    <td>Type 1-5</td>
    <td>--</td>
    <td>Not Found</td>
  </tr>
  <tr>
    <td>Reference Junction</td>
    <td>Type 1-5</td>
    <td>New Junction</td>
    <td>Precise or imprecise<br></td>
  </tr>
  <tr>
    <td>--</td>
    <td>--</td>
    <td>New Junction</td>
    <td>Novel</td>
  </tr>
</table>



## How to judge isoforms?

1.  **Ref Not found. This might happen in any of these cases:**
   a.  No exonic sequence overlap & no locus overlap
   b.  No exonic overlap but locus overlap (e.g. intronic containment
   c.  Minimal exonic overlap (but not “enough”)
   d.  Enough overlap but not the best match

2.  **The ref isoform was found with:**

	a.  Imprecise intronic chain
    b.  Precise intronic chain but wrong isoform boundaries
    c.  Precise intronic chain and precise isoform boundaries.
    d.  Same strand as reference
	e.  Opposite strand compared to the reference
	f.  unknown strand (unknown ref, unknown asm or unknown both)

> Note 2: Each category of the “found isoforms” can be further
> classified into
> 

	


3.  **Novel isoform**

    1.  Sharing insignificant sequence with reference isoform or sharing significant sequence but it is not the best one-to-one relationship.
    2.  Does not share any sequence with any reference isoform

> ***Note 1: The algorism of isoform matching:*** For each new assembled isoform, find all overlapping reference isoforms (overlapping means
> share exonic sequence). If it overlaps with many reference isoforms,
> assign this isoform to the reference isoform with longest overlap (as
> an absolute no of shared base pairs). If it shares the same length of
> sequence with multiple reference isoforms, choose the shortest
> reference isoform

<table>
  <tr>
    <th colspan="2">Reference assembly<br></th>
    <th colspan="3">New assembly<br></th>
  </tr>
  <tr>
    <th>Isoform id<br></th>
    <th>Gene</th>
    <th>Isoform id</th>
    <th>Gene</th>
    <th>Relationship</th>
  </tr>
  <tr>
    <td>Ref. isoform</td>
    <td>G1</td>
    <td>--</td>
    <td>--</td>
    <td>Not Found</td>
  </tr>
  <tr>
    <td>Ref. isoform</td>
    <td>G2</td>
    <td>Asm. Isoform</td>
    <td>NG1</td>
    <td>Imprecise, precise, or perfect<br></td>
  </tr>
  <tr>
    <td>--</td>
    <td>--</td>
    <td>Asm. Isoform</td>
    <td>NG2</td>
    <td>Novel-A</td>
  </tr>
  <tr>
    <td>Ref. isoform</td>
    <td>G3</td>
    <td>Asm. Isoform</td>
    <td>NG3</td>
    <td>Novel-B</td>
  </tr>
</table>

> **Note**: The same reference isoform can be assigned to many Asm. isoforms (many of them will have partial matches)
> 
> We can add extra step that test all the “Not Found” reference isoforms
> to define the exact subtype. Suggestion: We can make 2 versions of
> this table; one shows the “Novel” and “Not found” subtypes while the
> 2nd is concise and shows them as main category. This table will not
> have info about insignificant or not best matches so the last line in
> the above table will be removed.

## How to judge gene?

I do not think we need to test the gene. We can use the exon and transcript output tables to generate statistics about the genes

### Problems:

1.  **One of the problem that could face us while implementing the code is the fact of shared isoforms between different genes. If the some exon in isoform in gene x was shared with another isoform in gene y, an in the assembler, we found the isoform of gene x. Will we report that gene y isoform exists as well?**
    ***Solution***: The solution to this problem will be converting the problem from many-to-many relationship to one-to-one relationship. Instead of reporting all the isoforms overlapped between both reference and sequence, only isoform with the maximum overlapping is reported.
2.  **The second problem is the case where the best match for isoform in the reference, is not the same best match for the isoform in the assembly.**
    ***Solution***: A primary solution will be to start the comparison from one file and take the best match from this file’s side.

## Implementation:
To implement this code, we will need to convert the information in the GTF file to bunch of features for each judgment (exons, isoforms, and genes), and for each file as well (reference an assembly). By tracing the features from the reference in the assembly, we will get a hint on the sensitivity of the assembler, while tracing the features of the assembly will give us a hint on the assembler’s specificity.

## Pseudocode:

### Judging exons

> After converting the GTF for BED, each line will contain info about
> transcript (ID, no. of exons, start an end of the transcript, and
> start and end of each exon). However, the data about the genes is
> available only in the GTF file.

**We will use the GTF file to get the metadata about the genes (start and end), and to identify the gene where each transcript belong:**

- For each line in GTF where third column = gene

	- Add the gene id as a key in a genes\_dictionary, and set the value to start and end positions of this gene.

- For each line in GTF where third column = transcript

	-  Add the transcript as key in a transcripts\_dictionary and set the value to gene id, start, and end positions

	- we can also assign the transcript for its gene in genes dictionary. Will be easier when evaluating the genes.

- Now, to construct the table for exons in both ref. and asm., we will parse the lines from the BED file to get the exons of each transcript, and assign their metadata.

- For each line in BED file

	- Parse the start and end position for all exons

		- Identify transcript id

	- If first exon’s start position is the same as the gene’s start position

		- Exon’s type = 1

	- elif last exon’s end position is the same as gene’s end position

		- Exon’s type = 2

	- elif first exon’s start position is the same as transcript’s start position

		- Exon’s type = 4

	- elif first exon’s end position is the same as transcript’s end position

		- Exon’s type = 5

	- Else

		- Exon’s type = 3

	- Add the exon to exons dictionary and set the value to chr, start, end, transcript, gene and type

	- If exon already in exons dictionary, update the transcripts and genes lists, as well as exon’s type -if changed-

		- To make use of bedtools software on exonic data, we would import the relevant data from exons dictionary into a bed file, this includes chromosome, start, and end positions for both ref and asm.

		- intersectBed -a ref\_exon.bed -b asm\_exon.bed -wao &gt; output.bed

- wao asks bedtool to report all features with the amount of overlapping. Unfortunately, it comes as a number of overlapped base pairs rather than a fraction. We will then need to check for the overlap fraction to determine if an exon was found completely or partially.

- For each line in output.bed

	- If exon length == overlap

	- Overlap = 1

	- Add overlap fraction

- elif 0 &lt; exon length &lt; overlap

	- Overlap = overlap / exon length

	- Add overlap fraction

- Else \#in case of zero overlap

	- Add overlap fraction = Novel-ref

> Now, what was done is that for all asm. Exonic regions that overlap
> with the ref. Exonic regions have been reported. The next step is to
> check for asm. exons that weren’t reported, and assign their overlap
> to Novel-asm

- For each key in asm.exons dictionary

- If no overlap element is added yet OR overlap is less than a minimal overlap fraction

- Add overlap fraction = Novel-asm

### Judging exonic boundaries

> Having all the information about the exons reported in a table,
> judging boundaries will be a matter of exploring this table now.

For each line in exons\_output table

Define ref. And asm. Junctions positions and assign junction type based on exon type.

- If overlap == novel-ref

	- Junction\_relationshp = not found

	- If overlap == novel-asm

	- Junc\_rel = novel

- Else

	- If start and end position in both ref. And asm. Are equal

	- Junc\_rel = precise

- Else

	- Junc\_rel = ref\_junc\_start - asm\_junc\_start , ref\_junc\_end - asm\_junc\_end

### Judging isoforms

> - Judging the isoform will require some further manipulation, but i guess it will start by the same bedtools command, and the files are
> already prepared by converting the GTF to BED
> - intersectBed -a ref\_isoform.bed -b asm\_isoform.bed -wao &gt; isoform\_output.bed
> - The output of each line by this command will consist of the ref. Isoform id, the boundaries, no. of exons, exons boundaries, same data
> for the asm. Isoform, and finally the overlap.
> - The code for determining if an isoform is found or not can go as follows:

- For each line in isoform\_output.bed
	- If overlap = 0
		- Report not found
	- Elif overlap &lt; minimum number
		- Report novel asm. Isoform
	- Elif isoform start & end, exon numbers, and exon boundaries are the same in ref. And asm.
		-  OR exon numbers and exons boundaries are the same
		 - OR exon numbers are the same with imprecise exons boundaries
		- Report found
	- Else
		- In a dictionary for overlapped isoforms, store the ref. Isoform as key, and the asm. Isoform data
	- If the ref. Isoform already in the dict, append the new overlapped isoform data
		- For each key in overlapped isoform dictionary
			- Call the ref. Isoform data from the ref. Isoform dictionary
			- Identify the non-exonic parts \#the regions between the exons
			- If overlapped asm. Isoforms = 1
				- If all of its exons reside within the non-exonic region of the ref. Isoform
					- Report not found (or asm. novel)
				- Else
					- Report found
				- Else
					- Select the asm. Isoform with the longest overlap to be reported as found
					- And report the rest as asm. novel
				- For each overlapped asm. Isoform
					- If all its exons resides within the non-exonic region of the ref. Isoform
						- Report not found (or asm. novel)
					- Else
						- Assign scores for each feature in the ref. That’s found in asm.
						- Report the one with the highest score as found, and the rest as asm. Novel.
				- For each key in asm. Isoform
					- If isoform is not reported
						- Report asm. novel
