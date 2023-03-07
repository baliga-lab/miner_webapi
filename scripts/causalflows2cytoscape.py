#!/usr/bin/env python3
import os
import argparse
import MySQLdb
import json


DESCRIPTION = """causalflows2cytoscape.py - generate cytoscape data files"""
DBNAME = 'gbm_api'
MUTATION_REGULATOR_ROLES = { 1: 'down_regulates', 2: 'up_regulates'}
REGULATOR_REGULON_ROLES = { 1: 'activates', 2: 'represses' }

def dbconn():
    return MySQLdb.connect(host="localhost", user="admin", passwd="password",
                           db=DBNAME)


def make_regulator_files(conn, args):
    outpath = os.path.join(args.outdir, 'regulators')
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    with conn.cursor() as cur:
        cur.execute('select id,ensembl_id,preferred from genes where is_regulator=1')
        reg_counter = 1
        for reg_pk, reg_ensembl, reg_preferred in cur.fetchall():
            if reg_ensembl is not None:
                result = []
                # make json file for this regulator
                result.append({'data': {
                    #"id": "reg_%d" % reg_counter,
                    "id": reg_ensembl,
                    "name": reg_ensembl,
                    "classes": "tf"
                }})
                with conn.cursor() as cur2:
                    cur2.execute("""select r.name,regulon_regulator_spearman_r,m.name as mutation,cmft.name as cmftype,pw.name as pathway,mg.ensembl_id as mut_ensembl_id,mg.preferred as mut_preferred,tf.ensembl_id as tf_ensembl_id,tf.preferred as tf_preferred,regulon_mutation_regulator_role_id from cm_flows cmf join regulons r on cmf.regulon_id=r.id join mutations m on cmf.mutation_id=m.id join cm_flow_types cmft on cmf.cmf_type_id=cmft.id left outer join cmf_pathways as pw on pw.id=cmf.cmf_pathway_id left outer join genes mg on cmf.mutation_gene_id=mg.id join genes tf on cmf.regulator_id=tf.id where tf.ensembl_id=%s=%s""", [reg_ensembl, reg_ensembl])
                    mut_genes = set()
                    regulons = set()
                    regreg_edges = {}
                    mutreg_edges = {}
                    for row in cur2.fetchall():
                        (regulon, regulon_regulator_spearman_r, mutation, cmf_type, pathway,
                         mutgene_ensembl, mutgene_symbol,
                         regulator_ensembl, regulator_symbol,
                         mutation_regulator_role) = row
                        # 1. extract all mutation genes
                        mutgene = mutgene_symbol if mutgene_symbol is not None else mutgene_ensembl
                        mut_genes.add(mutgene)

                        # 2. extract all regulons
                        regulons.add(regulon)

                        # 3. extract the mutation-regulator edges
                        # 4. extract the regulator-regulon edges
                        regreg_edges['%s*%s' % (reg_ensembl, regulon)] = REGULATOR_REGULON_ROLES[1 if regulon_regulator_spearman_r < 0 else 2],
                        mutreg_edges['%s*%s' % (mutgene, reg_ensembl)] = MUTATION_REGULATOR_ROLES[mutation_regulator_role]
                    for mut in mut_genes:
                        result.append({'data': {
                            "id": mut,
                            "name": mut,
                            "classes": "mutation"
                        }})
                    for reg in regulons:
                        result.append({'data': {
                            "id": reg,
                            "name": reg,
                            "classes": "regulon"
                        }})
                    idcounter = 1
                    for key, role in regreg_edges.items():
                        reg_ensembl, regulon = key.split('*')
                        result.append({'data': {
                            "id": idcounter,
                            "source": reg_ensembl,
                            "target": regulon,
                            "classes": role
                        }})
                        idcounter += 1
                    for key, role in mutreg_edges.items():
                        mut, regulon = key.split('*')
                        result.append({'data': {
                            "id": idcounter,
                            "source": mut,
                            "target": regulon,
                            "classes": role
                        }})
                        idcounter += 1


                reg_counter += 1
                with open(os.path.join(outpath, "%s.json" % reg_ensembl), 'w') as outfile:
                    json.dump(result, outfile)
            else:
                raise Exception("no ENSEMBL ID !!!")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('outdir', help="output directory")
    args = parser.parse_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    conn = dbconn()
    try:
        make_regulator_files(conn, args)
    finally:
        conn.close()
