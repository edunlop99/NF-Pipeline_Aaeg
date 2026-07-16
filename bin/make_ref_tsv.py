import re, sys, os

if len(sys.argv) != 3:
    print("Usage: python make_ref_tsv.py <annotation.gtf_or_gff> <output.tsv>")
    sys.exit(1)

in_file = sys.argv[1]
out_file = sys.argv[2]

def parse_gtf_attrs(attr_string):
    gene_id      = re.search(r'gene_id "([^"]+)"', attr_string)
    gene_name    = re.search(r'gene_name "([^"]+)"', attr_string)
    gene_biotype = re.search(r'gene_biotype "([^"]+)"', attr_string)
    return (
        gene_id.group(1) if gene_id else None,
        gene_name.group(1) if gene_name else None,
        gene_biotype.group(1) if gene_biotype else None
    )

def parse_gff_attrs(attr_string):
    # handle both ID=gene:AAEL001 and ID=AAEL001 formats
    gene_id   = re.search(r'ID=(?:gene:)?([^;]+)', attr_string)
    gene_name = re.search(r'Name=([^;]+)', attr_string)
    biotype   = re.search(r'biotype=([^;]+)', attr_string)
    return (
        gene_id.group(1) if gene_id else None,
        gene_name.group(1) if gene_name else None,
        biotype.group(1) if biotype else None
    )

# gene-level feature types to look for
GTF_GENE_TYPES = {'gene'}
GFF_GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}

seen = set()
written = 0

with open(in_file) as f, open(out_file, 'w') as out:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue

        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        feature_type = fields[2]
        attr_string  = fields[8]

        # auto-detect format from attribute style
        is_gtf = '="' not in attr_string and '"' in attr_string

        if is_gtf:
            if feature_type not in GTF_GENE_TYPES:
                continue
            gid, gname, gbiotype = parse_gtf_attrs(attr_string)
        else:
            if feature_type not in GFF_GENE_TYPES:
                continue
            gid, gname, gbiotype = parse_gff_attrs(attr_string)

        if not gid or gid in seen:
            continue

        # fall back to gene_id if name/biotype missing
        gname    = gname    or gid
        gbiotype = gbiotype or 'unknown'

        seen.add(gid)
        out.write(f"{gid} {gname} {gbiotype}\n")
        written += 1

print(f"Done. Written {written} genes to {out_file}")