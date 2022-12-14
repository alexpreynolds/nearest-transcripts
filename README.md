# nearest-transcripts
Return nearest per-gene transcripts within N bases of some set of input regions ("sites")

## Setup

To install dependencies (`pyranges`, etc.):

```
$ pip install -r requirements.txt
```

## Usage

Example:

```
$ python find_nearest_transcripts.py > answer.txt
$ more answer.txt
Chromosome	Start	End	Transcript_Start	Transcript_End	Transcript_Strand	Transcript_ID	Transcript_Distance	Gene_ID
chr1	10000000	10000048	10032831	10180367	+	ENST00000253251.12	32784	ENSG00000130939.20
chr1	10000000	10000048	9981051	9985501	+	ENST00000496751.1	14500	ENSG00000173614.14
chr1	10000000	10000048	9930252	9943407	-	ENST00000377213.1	56594	ENSG00000162441.12
chr1	10000000	10000048	9850108	9910336	-	ENST00000400904.7	89665	ENSG00000178585.15
```

Edit `sites_local_fn` variable to point to a local file containing regions/sites of interest.