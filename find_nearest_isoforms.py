#!/usr/bin/env python

'''
find_nearest_isoforms.py

Return nearest per-gene isoforms within N bases of some set of regions
'''

import io
import os
import sys
import requests
import tempfile
import natsort as ns
import pandas as pd
import pyranges as pr

WINDOW_SIZE = 100_000

def main():
    '''
    set up temporary directory to store annotations and fake sites (TSS, TALE, etc.)
    '''
    temp_dir_handle = tempfile.TemporaryDirectory()
    temp_dir = temp_dir_handle.name

    '''
    grab annotations, e.g.:
    https://www.gencodegenes.org/human/release_42.html
    '''
    annotations_URI = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz'
    annotations_local_fn = os.path.join(temp_dir, 'gencode.v42.annotation.gtf.gz')
    if not os.path.exists(annotations_local_fn):
        sys.stderr.write('Getting annotations...\n')
        try:
            r = requests.get(annotations_URI)
            with open(annotations_local_fn, "wb") as ofh:
                b = io.BytesIO(r.content)
                ofh.write(b.getbuffer())
        except requests.exceptions.RequestException as e:
            raise SystemExit(e)
    
    '''
    read (for example) protein-coding transcripts into a Pandas dataframe (see: as_df=True option)
    https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/readers/index.html#pyranges.readers.read_gtf
    https://www.gencodegenes.org/pages/biotypes.html
    '''
    sys.stderr.write('Reading annotations into transcripts...\n')
    annotations_df = pr.read_gtf(annotations_local_fn, full=True, as_df=True, nrows=None, duplicate_attr=False)
    transcripts_df = (annotations_df
        .query("(gene_type == 'protein_coding')")
        .query("(Feature == 'transcript')")
        .set_index('gene_id')
    )

    '''
    convert Pandas dataframe into PyRanges object
    '''
    transcripts = pr.PyRanges(transcripts_df.reset_index())

    '''
    generate file of fake sites (optional)
    '''
    sys.stderr.write('Generating fake sites...\n')
    sites_local_fn = os.path.join(temp_dir, 'sites.bed')
    if not os.path.exists(sites_local_fn):
        fake_site = 'chr1\t10000000\t10000048\n' # some random 48nt region
        with open(sites_local_fn, 'w') as ofh:
            ofh.write(fake_site)

    '''
    read sites into PyRanges object
    NOTE: make it a global variable, so that it can be used later on in per_gene_nearest_transcript_filter()
    https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/index.html#pyranges.read_bed
    '''
    sys.stderr.write('Reading sites into PyRanges object...\n')
    global sites
    sites = pr.read_bed(sites_local_fn, as_df=False, nrows=None)

    '''
    1) extend boundaries of sites by WINDOW_SIZE constant (e.g., 100kb)
    2) extended regions might overlap, and are merged, ignoring strandedness (if present)
    '''
    sys.stderr.write('Extending site regions by {}nt...\n'.format(WINDOW_SIZE))
    sites_extended = sites.extend(WINDOW_SIZE).merge(strand=False)

    '''
    find all transcripts that overlap extended sites
    NOTE: ignores strandedness of sites; if needed, adjust sites file with strand data, refer
          to the API for options to deal with strandedness, and test the results that come back
    https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/index.html#pyranges.PyRanges.join
    '''
    sys.stderr.write('Finding transcripts that overlap extended site regions...\n')
    transcripts_overlapping_extended_sites = transcripts.join(sites_extended, apply_strand_suffix=False)

    '''
    filter each per-gene set of transcripts gene for those nearest to each site, and
    convert it to a Pandas dataframe of results
    '''
    sys.stderr.write('Finding per-gene transcripts nearest to original sites...\n')
    results = transcripts_overlapping_extended_sites.apply(per_gene_nearest_transcript_filter).df

    '''
    write sorted, formatted results to standard output
    '''
    sys.stderr.write('Writing results to stdout...\n')
    results = results.rename(columns={
        "Start_b": "Transcript_Start",
        "End_b": "Transcript_End",
        "Strand": "Transcript_Strand",
        "transcript_id": "Transcript_ID",
        "Distance": "Transcript_Distance",
        "gene_id": "Gene_ID",
    })
    results = results.sort_values(by=["Chromosome", "Start", "End"])
    results = results.reindex(
        index=ns.order_by_index(
            results.index, 
            ns.index_natsorted(zip(
                results.Chromosome, 
                results.Start, 
                results.End
            ))
        )
    )
    results = results.reset_index(drop=True)
    output = io.StringIO()
    results.to_csv(output, sep="\t", index=False, header=True)
    sys.stdout.write(f"{output.getvalue()}")

    '''
    cleanup
    '''
    sys.stderr.write('Cleaning up temporary directory...\n')
    temp_dir_handle.cleanup()

'''
1) group overlapping transcripts by each gene_id value
2) find the nearest transcript to each site
3) return a PyRanges object containing the site and the ID of the nearest gene and transcript pair
NOTE: using overlap=False to filter out transcripts that overlap the site
https://pyranges.readthedocs.io/en/latest/autoapi/pyranges/index.html#pyranges.PyRanges.nearest
'''
def per_gene_nearest_transcript_filter(df):
    answer = pd.DataFrame()
    grouped_per_gene_id_df = df.groupby(by='gene_id')
    for per_gene_id, per_gene_id_group in grouped_per_gene_id_df:
        grouped_per_gene_id_pr = pr.PyRanges(per_gene_id_group)
        grouped_per_gene_id_column_subset = grouped_per_gene_id_pr.drop(drop=None)
        setattr(grouped_per_gene_id_column_subset, 'transcript_id', per_gene_id_group['transcript_id'])
        nearest_per_gene_transcript_to_each_site = sites.nearest(grouped_per_gene_id_column_subset, strandedness=None, overlap=False, how=None, nb_cpu=1) #.drop(drop=None)
        setattr(nearest_per_gene_transcript_to_each_site, 'gene_id', per_gene_id)
        if not answer.empty:
            answer = pd.concat([answer, nearest_per_gene_transcript_to_each_site.df], ignore_index=True)
        else:
            answer = nearest_per_gene_transcript_to_each_site.df
    return answer

if __name__ == "__main__":
    main()