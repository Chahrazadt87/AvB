"""Supplementary table of cinquedea (hydrolase 2090) homologs.

Sheets:
    README                 column descriptions
    Sequence homologs      MMseqs2 search against GTDB r214 and UniProtKB 2023_05 (Figure 3)
    Structural homologs    Foldseek search against AFDB50 (Foldseek figures)

Run from the repo root:
    python src/homology_search/make_supplementary_table.py
"""
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

sys.path.insert(0, 'src/foldseek_2090')
from foldseek_2090 import W as FOLDSEEK_DIR, load_pfam, load_search_hits

DATA = Path('data')
SEARCH = DATA / 'hydrolase_search'
GTDB = Path('/Users/rs1521/workspace/data/gtdb_r214.1')
OUT = DATA / 'supplementary'
QUERY_ID = 'pgaptmp_002090_1'

CONTIG_CACHE = SEARCH / 'gtdb_contig_to_assembly.json'


# --- Genome accessions of GTDB proteins -------------------------------------------------------

def gtdb_genomes():
    m = pd.concat([
        pd.read_csv(GTDB / f, sep='\t', low_memory=False,
                    usecols=['accession', 'ncbi_genbank_assembly_accession', 'ncbi_wgs_master', 'gtdb_taxonomy'])
        for f in ['ar53_metadata_r214.tsv', 'bac120_metadata_r214.tsv']
    ])
    m['number'] = m.ncbi_genbank_assembly_accession.str.extract(r'GC[AF]_(\d+)')[0]
    m['wgs'] = m.ncbi_wgs_master.str.extract(r'^([A-Z]{4,6})\d')[0]
    return m


def contig_to_assemblies(contigs):
    """NCBI Datasets lookup of the assemblies containing each contig (cached).

    Only resolves complete sequences (e.g. chromosomes), not WGS contigs.
    """
    cache = json.loads(CONTIG_CACHE.read_text()) if CONTIG_CACHE.exists() else {}
    cache = {k: v for k, v in cache.items() if isinstance(v, list)}
    for c in sorted(set(contigs) - set(cache)):
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/sequence_accession/{c}/sequence_assemblies'
        r = requests.get(url, timeout=60)
        cache[c] = r.json().get('accessions', []) if r.ok else []
        time.sleep(0.35)
    CONTIG_CACHE.write_text(json.dumps(cache, indent=1))
    return cache


def genome_accession(row, contig_assemblies, gtdb):
    """Assembly accession of the genome encoding a protein."""
    for col in ['id', 'db_proka_id']:
        v = row[col]
        if isinstance(v, str) and '@' in v:
            return v.split('@')[1]
    if isinstance(row['gtdb_id'], str):
        contig = row['gtdb_id'].rsplit('_', 1)[0]
        species = gtdb.gtdb_taxonomy.str.endswith(f's__{row["gtdb_species"]}')
        # WGS contigs: the accession prefix identifies the WGS project of the genome.
        wgs = re.match(r'^(?:NZ_)?([A-Z]{4,6})\d', contig)
        if wgs:
            hits = gtdb[(gtdb.wgs == wgs.group(1)) & species]
            if len(hits):
                return re.sub(r'^(RS|GB)_', '', hits.iloc[0].accession)
        candidates = contig_assemblies.get(contig, [])
        numbers = {re.search(r'GC[AF]_(\d+)', a).group(1): a for a in candidates}
        # Prefer the assembly that is part of GTDB r214, matching on species.
        hits = gtdb[gtdb.number.isin(numbers) & species]
        if len(hits):
            return re.sub(r'^(RS|GB)_', '', hits.iloc[0].accession)
        if candidates:
            return sorted(candidates)[0]
        # Last resort: the species is represented by a single genome in GTDB r214.
        if species.sum() == 1:
            return re.sub(r'^(RS|GB)_', '', gtdb[species].iloc[0].accession)
    return None


def ncbi_protein_accession(row):
    for col in ['id', 'db_proka_id']:
        v = row[col]
        if isinstance(v, str) and '@' in v and not re.match(r'.+\.\d+_\d+$', v.split('@')[0]):
            return v.split('@')[0]
    return None


# --- Pfam ---------------------------------------------------------------------------------------

def format_pfam(df, id_col, name_col, acc_col, start_col, end_col, evalue_col):
    df = df.sort_values([id_col, start_col])
    return df.groupby(id_col).apply(lambda g: '; '.join(
        f'{r[name_col]} ({r[acc_col].split(".")[0]}; {r[start_col]}-{r[end_col]}; E={r[evalue_col]:.1e})'
        for _, r in g.iterrows()
    ))


# --- Sheets -------------------------------------------------------------------------------------

def sequence_homologs():
    d = pd.read_csv(SEARCH / 'search_output.csv')
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(SEARCH / 'search_output.fasta', 'fasta')}
    tree_ids = {r.id for r in SeqIO.parse(DATA / 'hydrolase_tree' / 'alignment_final.fasta', 'fasta')}

    gtdb = gtdb_genomes()
    need = d[d.gtdb_id.notna() & ~d.id.str.contains('@') & ~d.db_proka_id.fillna('').str.contains('@')]
    contig_assemblies = contig_to_assemblies(need.gtdb_id.str.rsplit('_', n=1).str[0])

    pfam = pd.read_csv(SEARCH / 'search_output.pfam.csv')
    pfam_str = format_pfam(pfam, 'protein_id', 'hmm_query', 'hmm_accession', 'start', 'end', 'evalue')

    t = pd.DataFrame({
        'protein_id': d.id,
        'source_database': [
            'GTDB r214 + UniProtKB' if isinstance(g, str) and isinstance(u, str)
            else 'GTDB r214' if isinstance(g, str) else 'UniProtKB'
            for g, u in zip(d.gtdb_id, d.uniprot_id)
        ],
        'uniprot_accession': d.uniprot_id,
        'ncbi_protein_accession': d.apply(ncbi_protein_accession, axis=1),
        'gtdb_gene_id': d.gtdb_id,
        'genome_accession': d.apply(genome_accession, axis=1, contig_assemblies=contig_assemblies, gtdb=gtdb),
        'domain': d.domain,
        **{f'gtdb_{r}': d[f'gtdb_{r}'] for r in ['phylum', 'class', 'order', 'family', 'genus', 'species']},
        **{f'ncbi_{r}': d[f'ncbi_{r}'] for r in ['phylum', 'class', 'order', 'family', 'genus', 'species']},
        'mmseqs_evalue': d.evalue,
        'mmseqs_bitscore': d.bits,
        'target_start': d.tstart,
        'target_end': d.tend,
        'pfam_domains': d.id.map(pfam_str),
        'in_protein_tree': d.id.isin(tree_ids).map({True: 'yes', False: 'no (>50% gaps after trimming)'}),
        'sequence': d.id.map(seqs),
    })

    query = {
        'protein_id': QUERY_ID, 'source_database': 'This study (query)', 'domain': 'Archaea',
        'gtdb_species': 'Haloferax larsenii s5a-1', 'in_protein_tree': 'yes',
        'pfam_domains': pfam_str.get(QUERY_ID), 'sequence': seqs[QUERY_ID],
    }
    return pd.concat([pd.DataFrame([query]), t], ignore_index=True)


def structural_homologs(sequence_table):
    hits = load_search_hits()
    table = pd.read_csv(f'{FOLDSEEK_DIR}/hits_annotated.tsv', sep='\t', index_col='id')
    no_euk = pd.read_csv(f'{FOLDSEEK_DIR}/hits_annotated_no_eukaryotes.tsv', sep='\t', index_col='id')
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(f'{FOLDSEEK_DIR}/all_hits.fasta', 'fasta')}

    pfam = load_pfam()
    pfam_str = format_pfam(pfam, 'target', 'query', 'accession', 'env_from', 'env_to', 'ievalue')

    uniprot = hits.id.str.extract(r'AF-(.+)-F1')[0]
    t = pd.DataFrame({
        'afdb_id': hits.id,
        'uniprot_accession': uniprot,
        'afdb_v6_model_available': hits.id.isin(table.index[table.has_structure]).map({True: 'yes', False: 'no (obsolete UniProtKB entry)'}),
        'ncbi_taxid': hits.taxid,
        'organism': hits.taxname,
        'domain': hits.domain.fillna('Metagenome (unassigned)'),
        'ncbi_phylum': hits.phylum,
        'ncbi_class': hits['class'],
        'foldseek_rank': hits['rank'],
        'foldseek_evalue': hits.evalue,
        'foldseek_bitscore': hits.bits,
        'foldseek_qtmscore': hits.qtmscore,
        'foldseek_sequence_identity': hits.fident,
        'foldseek_alignment_length': hits.alnlen,
        'tm_score_to_cinquedea': hits.id.map(table.tm_to_2090),
        'pcoa_1': hits.id.map(table.PC1),
        'pcoa_2': hits.id.map(table.PC2),
        'pcoa_1_no_eukaryotes': hits.id.map(no_euk.PC1),
        'pcoa_2_no_eukaryotes': hits.id.map(no_euk.PC2),
        'pfam_domains': hits.id.map(pfam_str),
        'pfam_category': hits.id.map(table.category),
        'in_sequence_homologs': uniprot.isin(set(sequence_table.uniprot_accession.dropna())).map({True: 'yes', False: 'no'}),
        'sequence': hits.id.map(seqs),
    })
    return t


README = [
    ('Sequence homologs', None),
    ('protein_id', 'Identifier used in the homology search and protein tree (Figure 3A)'),
    ('source_database', 'Database(s) in which the homolog was found (MMseqs2, -s 7.5): GTDB release 214 and/or UniProtKB release 2023_05; GTDB hits were deduplicated against UniProtKB (see Methods)'),
    ('uniprot_accession', 'UniProtKB accession'),
    ('ncbi_protein_accession', 'NCBI RefSeq/GenBank protein accession, where available'),
    ('gtdb_gene_id', 'GTDB r214 gene identifier (nucleotide accession of the contig followed by the gene number)'),
    ('genome_accession', 'NCBI assembly accession of the genome encoding the protein (GTDB hits only)'),
    ('domain', 'Domain of life'),
    ('gtdb_phylum ... gtdb_species', 'GTDB r214 taxonomy (GTDB hits)'),
    ('ncbi_phylum ... ncbi_species', 'NCBI taxonomy (UniProtKB hits)'),
    ('mmseqs_evalue, mmseqs_bitscore', 'MMseqs2 E-value and bit score of the alignment with cinquedea'),
    ('target_start, target_end', 'Aligned region of the homolog'),
    ('pfam_domains', 'Pfam (release 37) domains predicted with HMMER hmmsearch using gathering thresholds (--cut_ga): name (accession; envelope start-end; independent E-value)'),
    ('in_protein_tree', 'Whether the protein is included in the tree of Figure 3A'),
    ('sequence', 'Amino acid sequence'),
    ('', None),
    ('Structural homologs', None),
    ('afdb_id', 'AlphaFold Protein Structure Database identifier'),
    ('uniprot_accession', 'UniProtKB accession'),
    ('afdb_v6_model_available', 'Whether an AFDB v6 model could be retrieved for all-vs-all comparisons (Foldseek search was performed against AFDB50 v4)'),
    ('ncbi_taxid, organism, domain, ncbi_phylum, ncbi_class', 'NCBI taxonomy reported by Foldseek'),
    ('foldseek_*', 'Foldseek search of the cinquedea AlphaFold3 model against AFDB50 (--cluster-search 1): rank (by bit score), E-value, bit score, TM-score normalised by query length, sequence identity and alignment length'),
    ('tm_score_to_cinquedea', 'Structural similarity to cinquedea from the all-vs-all Foldseek comparison: mean of TM-scores normalised by each protein length, averaged over both search directions'),
    ('pcoa_1, pcoa_2', 'Coordinates in the principal coordinates analysis of all structural homologs'),
    ('pcoa_1_no_eukaryotes, pcoa_2_no_eukaryotes', 'Coordinates in the principal coordinates analysis excluding eukaryotic homologs'),
    ('pfam_domains', 'Pfam (release 38) domains predicted with HMMER hmmsearch using gathering thresholds (--cut_ga): name (accession; envelope start-end; independent E-value)'),
    ('pfam_category', 'Pfam category used in the figures'),
    ('in_sequence_homologs', 'Whether the protein is also among the sequence homologs (UniProtKB accession match)'),
    ('sequence', 'Amino acid sequence (last archived UniProtKB version for obsolete entries)'),
]


if __name__ == '__main__':
    OUT.mkdir(exist_ok=True)
    seq = sequence_homologs()
    struct = structural_homologs(seq)

    readme = pd.DataFrame(README, columns=['column', 'description'])
    path = OUT / 'Table_S4_cinquedea_homologs.xlsx'
    with pd.ExcelWriter(path) as w:
        readme.to_excel(w, sheet_name='README', index=False)
        seq.to_excel(w, sheet_name='Sequence homologs', index=False)
        struct.to_excel(w, sheet_name='Structural homologs', index=False)
    seq.to_csv(OUT / 'Table_S4_sequence_homologs.tsv', sep='\t', index=False)
    struct.to_csv(OUT / 'Table_S4_structural_homologs.tsv', sep='\t', index=False)
    print(path, len(seq), len(struct))
