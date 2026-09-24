"""Helpers for the Foldseek analysis of hydrolase 2090 (cinquedea) hits."""
import numpy as np
import pandas as pd

W = 'data/foldseek_2090'
QUERY = 'hydrolase_2090_AF3'

DOMAIN_COLORS = {
    'Archaea': '#f58231',
    'Bacteria': '#4363d8',
    'Eukaryota': '#a9a9a9',
}

PFAM_ORDER = ['DUF900', 'DUF726', 'DUF676', 'Other α/β hydrolase (CL0028)',
              'Periplasmic binding protein (CL0144)', 'P-loop NTPase (CL0023)', 'Other Pfam', 'No Pfam hit']
PFAM_COLORS = {
    'DUF900': '#e6194b',
    'DUF726': '#911eb4',
    'DUF676': '#f032e6',
    'Other α/β hydrolase (CL0028)': '#3cb44b',
    'Periplasmic binding protein (CL0144)': '#808000',
    'P-loop NTPase (CL0023)': '#000075',
    'Other Pfam': '#9a6324',
    'No Pfam hit': '#d9d9d9',
}


def load_search_hits(path='data/hydrolase_2090_AF3_foldseek_results.tsv'):
    """Original Foldseek search of 2090 against AFDB50 (October 2024)."""
    d = pd.read_csv(path, sep='\t')
    d['id'] = d.target.str.replace('-model_v4', '', regex=False)
    for rank in ['d', 'p', 'c', 'o', 'f', 'g']:
        d[rank] = d.taxlineage.str.extract(rf'(?:^|;){rank}_([^;]+)')
    d = d.rename(columns={'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus'})
    d['rank'] = np.arange(1, len(d) + 1)
    return d


def load_pfam(path=f'{W}/all_hits.pfam.domtblout.txt',
              clans='/Users/rs1521/workspace/data/Pfam_38/Pfam-A.clans.tsv'):
    cols = ['target', 'query', 'accession', 'ievalue', 'env_from', 'env_to']
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            x = line.split()
            rows.append((x[0], x[3], x[4].split('.')[0], float(x[12]), int(x[19]), int(x[20])))
    pfam = pd.DataFrame(rows, columns=cols)
    clan = pd.read_csv(clans, sep='\t', header=None, names=['accession', 'clan', 'clan_name', 'name', 'description'])
    return pfam.merge(clan[['accession', 'clan']], on='accession', how='left')


def pfam_category(names, clans):
    names, clans = set(names), set(clans)
    for n in ['DUF900', 'DUF726', 'DUF676']:
        if n in names:
            return n
    if 'CL0028' in clans:
        return 'Other α/β hydrolase (CL0028)'
    if 'CL0144' in clans:
        return 'Periplasmic binding protein (CL0144)'
    if 'CL0023' in clans:
        return 'P-loop NTPase (CL0023)'
    return 'Other Pfam' if names else 'No Pfam hit'


def annotate_hits(hits, pfam):
    per_seq = pd.DataFrame({
        target: {
            'pfam_domains': ';'.join(dict.fromkeys(g['query'])),
            'category': pfam_category(g['query'], g['clan'].dropna()),
        }
        for target, g in pfam.sort_values(['target', 'env_from']).groupby('target')
    }).T
    out = hits.merge(per_seq, left_on='id', right_index=True, how='left')
    out['pfam_domains'] = out.pfam_domains.fillna('')
    out['category'] = out.category.fillna('No Pfam hit')
    return out


def load_all_vs_all(path=f'{W}/all_vs_all.tsv'):
    return pd.read_csv(
        path, sep='\t', header=None,
        names=['query', 'target', 'fident', 'alnlen', 'qlen', 'tlen', 'evalue', 'bits',
               'qtmscore', 'ttmscore', 'alntmscore'],
    )


def similarity_matrix(ava, ids, score='mean_tm'):
    """Symmetric structural similarity; pairs without an alignment get 0.

    mean_tm: average of the TM-scores normalised by each chain's length,
    then averaged over both search directions.
    """
    idx = {k: i for i, k in enumerate(ids)}
    ava = ava[ava['query'].isin(idx) & ava['target'].isin(idx)]
    if score == 'mean_tm':
        s = (ava.qtmscore + ava.ttmscore) / 2
    else:
        s = ava[score]
    n = len(ids)
    m = np.zeros((n, n))
    i = ava['query'].map(idx).to_numpy()
    j = ava['target'].map(idx).to_numpy()
    m[i, j] = s.to_numpy()
    m = (m + m.T) / 2
    np.fill_diagonal(m, 1)
    return m


def pcoa(d, k=3):
    """Classical multidimensional scaling of a distance matrix."""
    n = len(d)
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    coords = vecs[:, :k] * np.sqrt(np.maximum(vals[:k], 0))
    explained = vals[:k] / vals[vals > 0].sum()
    return coords, explained
