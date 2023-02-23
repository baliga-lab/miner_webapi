#!/usr/bin/env python3

import os
import argparse
import MySQLdb
import pandas as pd
import numpy as np
import json
from collections import defaultdict


DESCRIPTION = "import_gbm - import filtered causal results to GBM backend"

"""
Tables in mm_api_v2:


ignored for now (because patients are needed):

patient_tf
patients
bicluster_patients
bicluster_boxplot_data

bc_tf_roles: 1 = activates, 2 = represses
bc_mutation_tf_roles: 1 = down-regulates 2 = up-regulates
"""

def fill_roles(conn):
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_tf_roles')
        # only insert if not exist
        if cur.fetchone()[0] == 0:
            cur.execute("insert into bc_tf_roles (id, name) values (1, 'activates')")
            cur.execute("insert into bc_tf_roles (id, name) values (2, 'represses')")
            cur.execute("insert into bc_mutation_tf_roles (id, name) values (1, 'down-regulates')")
            cur.execute("insert into bc_mutation_tf_roles (id, name) values (2, 'up-regulates')")
            conn.commit()

def import_tfs(conn, df):
    """tfs: id: int, name: string, cox_hazard_ratio: float
    In the import document this is the Regulator"""
    tfs = set()
    for index, row in df.iterrows():
        tfs.add(row['Regulator'])
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from tfs')
        num_tfs = cur.fetchone()[0]
        if num_tfs > 0:
            print('TFs found, reading from database')
            cur.execute('select id,name from tfs')
            for pk, name in cur.fetchall():
                result[name] = pk
        else:
            print('TFs not found, inserting into database')
            for tf in tfs:
                cur.execute('insert into tfs (name) values (%s)', [tf])
                result[tf] = cur.lastrowid
            conn.commit()
    return result


def import_mutations(conn, df):
    """mutations: id: int, name: string"""
    mutations = sorted(set(df["Mutation"]))

    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from mutations')
        if cur.fetchone()[0] > 0:
            print("Mutations found, reading from database")
            cur.execute('select id,name from mutations')
            for pk, name in cur.fetchall():
                result[name] = pk
        else:
            print("Mutations not found, inserting into database")
            for name in mutations:
                cur.execute('insert into mutations (name) values (%s)', [name])
                result[name] = cur.lastrowid
            conn.commit()
    return result


def import_genes(conn, genes, ens2pref, ens2entrez):
    """genes: id: int, ensembl_id: string, entrez_id: string, preferred: string
    In the import document this is the Regulator
    """
    print("# unique genes found: %d" % len(genes))

    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from genes')
        if cur.fetchone()[0] > 0:
            print("Genes found reading from database")
            cur.execute('select id,ensembl_id,entrez_id,preferred from genes')
            for pk, ens, entrez, pref in cur.fetchall():
                result[ens] = (pk, entrez, pref)
        else:
            print("Genes not found, inserting into database")
            for ens in sorted(genes):
                try:
                    entrez = ens2entrez[ens]
                except:
                    entrez = None
                try:
                    pref = ens2pref[ens]
                except:
                    pref = None
                cur.execute('insert into genes (ensembl_id,entrez_id,preferred) values (%s,%s,%s)',
                            [ens, entrez, pref])
                result[ens] = (cur.lastrowid, entrez, pref)
            conn.commit()
    return result


def import_regulons(conn, regulon_map, cox_map, mutations):
    """biclusters: id, name, cox_hazard_ratio"""
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from biclusters')
        if cur.fetchone()[0] > 0:
            print('Regulons found, reading from database')
            cur.execute('select b.id,b.name,b.cox_hazard_ratio from biclusters b')
            for pk, regulon, hr in cur.fetchall():
                result[regulon] = (pk, hr)
        else:
            print('Regulons not found, inserting into database')
            for regulon in regulon_map.keys():
                cox_hr = cox_map[regulon]['HazardRatio']
                cur.execute('select count(*) from biclusters where name=%s', [regulon])
                if cur.fetchone()[0] == 0:
                    cur.execute('insert into biclusters (name,cox_hazard_ratio) values (%s,%s)',
                                [regulon, cox_hr])
                    result[regulon] = (cur.lastrowid, cox_hr)
            conn.commit()
    return result


def import_mutation_regulator(conn, df, regulons, mutations, tfs):
    """table: bc_mutation_tf, field: MutationRegulatorEdge
    id, bicluster_id, mutation_id, tf_id, role"""
    """bc_tf_roles: 1 = activates, 2 = represses
    bc_mutation_tf_roles: 1 = down-regulates 2 = up-regulates"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_mutation_tf')
        if cur.fetchone()[0] == 0:
            print('no mutation_regulator edges found, insert into database')

            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulator = row['Regulator']
                mutation = row['Mutation']
                edge = row['MutationRegulatorEdge']
                regulon_id, cox_hr1  = regulons[regulon]
                regulator_id = tfs[regulator]
                mutation_id = mutations[mutation]
                if edge < 0:
                    role = 1  # down-regulates
                else:
                    role = 2  # up-regulates
                cur.execute('insert into bc_mutation_tf (bicluster_id,mutation_id,tf_id,role) values (%s,%s,%s,%s)',
                            [regulon_id, mutation_id, regulator_id, role])
            conn.commit()
        else:
            print('skip inserting mutation regulator edges')

"""
Index(['Unnamed: 0', 'Mutation', 'Regulator', 'Regulon',
       'MutationRegulatorEdge', '-log10(p)_MutationRegulatorEdge',
       'RegulatorRegulon_Spearman_R', 'RegulatorRegulon_Spearman_p-value',
       'Regulon_stratification_t-statistic',
       '-log10(p)_Regulon_stratification',
       'Fraction_of_edges_correctly_aligned', 'TranscriptionalProgram',
       'HazardRatio', 'HazardRatioPval', 'RegulonGenes', 'DrugEnrichment',
       'HallmarksEnrichment', 'LinHallmarks', 'TargetClassEnrichment',
       'TargetClass_p-value', 'MechanismOfActionEnrichment',
       'MechanismOfAction_p-value', 'mirv_miRNA', 'PITA_miRNA', 'PITA_pval',
       'TargetScan_miRNA', 'TargetScan_pval'],
      dtype='object')

"""

def import_regulon_regulator(conn, df, regulons, tfs):
    """table: bc_tf (bicluster_id,tf_id,role)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_tf')
        if cur.fetchone()[0] == 0:
            print('insert regulon-regulator into database')
            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulator = row['Regulator']
                regulon_id = regulons[regulon][0]
                regulator_id = tfs[regulator]
                spearman_r = row['RegulatorRegulon_Spearman_R']
                if spearman_r > 0:
                    role = 1  # 1 activates
                else:
                    role = 2  # 2 represses
                #print(row['RegulatorRegulon_Spearman_p-value'])
                cur.execute('select count(*) from bc_tf where bicluster_id=%s and tf_id=%s',
                            [regulon_id, regulator_id])
                if cur.fetchone()[0] == 0:
                    # only insert once
                    cur.execute('insert into bc_tf (bicluster_id,tf_id,role) values (%s,%s,%s)',
                                [regulon_id, regulator_id, role])
            conn.commit()
        else:
            print('regulon regulator relations exist, skip')


def import_regulon_genes(conn, df, regulon_map, regulons, genes):
    """table: bicluster_genes (bicluster_id, gene_id)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from bicluster_genes')
        if cur.fetchone()[0] == 0:
            print('insert regulon genes into database')
            for regulon, regulon_genes in regulon_map.items():
                regulon_id = regulons[regulon][0]
                for gene in regulon_genes:
                    gene_id = genes[gene][0]
                    cur.execute('select count(*) from bicluster_genes where bicluster_id=%s and gene_id=%s',
                                [regulon_id, gene_id])
                    if cur.fetchone()[0] == 0:  # prevent double insertions in the same transaction
                        cur.execute('insert into bicluster_genes(bicluster_id,gene_id) values (%s,%s)',
                                    [regulon_id, gene_id])
            conn.commit()
        else:
            print('regulon genes exist, skip')


def import_transcriptional_programs(conn, regulons, args):
    with open(os.path.join(args.indir, 'transcriptional_programs.json')) as infile:
        regulon2programs = json.load(infile)
    print('import transcriptional programs')
    result = {}
    programs = {}
    with conn.cursor() as cur:
        for regulon_number, reg_programs in regulon2programs.items():
            # add to database if necessary
            regulon = 'R-%s' % regulon_number
            regulon_id = regulons[regulon][0]
            for program in reg_programs:
                progname = "P-%s" % program
                if progname not in programs:
                    cur.execute('insert into trans_programs (name) values (%s)', [progname])
                    programs[progname] = cur.lastrowid
                cur.execute('insert into bicluster_programs (bicluster_id, program_id) values (%s,%s)', [regulon_id, programs[progname]])
        conn.commit()
    return programs

def import_regulons_programs_genes_disease_mapping(conn,
                                                   regulons, programs, genes,
                                                   args):
    df = pd.read_csv(os.path.join(args.indir, 'regulons_programs_genes_disease_mapping.csv'),
                     sep=',', header=0)
    with conn.cursor() as cur:
        cur.execute('select count(*) from bc_program_genes')
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            return  # Skip, we already imported
        for index, row in df.iterrows():
            regulon = 'R-%d' % row['Regulon_ID']
            gene = row['Gene']
            program = 'P-%d' % row['Programs']
            regulon_id = regulons[regulon][0]
            gene_id = genes[gene][0]
            program_id = programs[program]
            is_disease_relevant = row['is_disease_relevant']
            print("GENEID: ", gene_id)
            print("REGID: ", regulon_id)
            print("PROGID: ", program_id)
            cur.execute('insert into bc_program_genes (bicluster_id,program_id,gene_id,is_disease_relevant) values (%s,%s,%s,%s)', [regulon_id, program_id, gene_id, 1 if is_disease_relevant != 'False' else 0])
        conn.commit()


def import_type_table(conn, table, values):
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from %s' % table)
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            cur.execute('select id, name from %s' % table)
            for pk, name in cur.fetchall():
                result[name] = pk
        else:
            for value in values:
                cur.execute('insert into ' + table + ' (name) values (%s)', [value])
                result[value] = cur.lastrowid
            conn.commit()

    return result

def import_drug_types(conn, values):
    return import_type_table(conn, 'drug_types', values)


def import_action_types(conn, values):
    return import_type_table(conn, 'action_types', values)


def import_mechanisms_of_action(conn, values):
    return import_type_table(conn, 'mechanisms_of_action', values)


def import_target_type_models(conn, values):
    return import_type_table(conn, 'target_type_models', values)


def import_drugs(conn, df, drug_types, action_types, mechanisms_of_action, target_type_models):
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from drugs')
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            cur.execute('select id,name from drugs')
            for pk, name in cur.fetchall():
                result[name] = pk
        else:
            for index, row in df.iterrows():
                drug_name = row['Drug']
                if drug_name in result:
                    continue  # already inserted, skip
                approved_symbol = row['approvedSymbol']
                is_approved = row['isApproved'] == 'TRUE'
                max_trial_phase = row['maximumClinicalTrialPhase']
                if np.isnan(max_trial_phase):
                    max_trial_phase = None
                max_phase_gbm = row['maxPhaseGBM']
                if np.isnan(max_phase_gbm):
                    max_phase_gbm = None
                drug_type_id = drug_types[row['drugType']]
                action_type_id = action_types[row['actionType']]
                mechanism_of_action_id = mechanisms_of_action[row['mechanismOfAction']]
                target_type_model_id = target_type_models[row['TargetTypeModel']]
                cur.execute('insert into drugs (name,approved_symbol,approved,max_trial_phase,max_phase_gbm,drug_type_id,action_type_id,target_type_model_id,mechanism_of_action_id) values (%s,%s,%s,%s,%s,%s,%s,%s,%s)', [drug_name, approved_symbol, 1 if is_approved else 0, max_trial_phase, max_phase_gbm, drug_type_id, action_type_id, mechanism_of_action_id, target_type_model_id])
                result[drug_name] = cur.lastrowid
            conn.commit()

    return result

def import_drug_targets(conn, df, drugs, genes):
    with conn.cursor() as cur:
        cur.execute('select count(*) from drug_targets')
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            return  # skip
        seen = set()
        for index, row in df.iterrows():
            drug_id = drugs[row['Drug']]
            try:
                target_id = genes[row['targets']][0]
                key = '%d-%d' % (drug_id, target_id)
                if key not in seen:
                    cur.execute('insert into drug_targets (drug_id, gene_id) values (%s,%s)', [drug_id, target_id])
                    seen.add(key)
            except KeyError:
                pass  ## skip
        conn.commit()


def import_drug_regulons(conn, df, drugs, regulons):
    with conn.cursor() as cur:
        cur.execute('select count(*) from drug_regulons')
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            return  # skip
        seen = set()
        for index, row in df.iterrows():
            drug_id = drugs[row['Drug']]
            regulon_id = regulons['R-%d' % row['Regulon_ID']][0]
            key = '%d-%d' % (drug_id, regulon_id)
            if key not in seen:
                cur.execute('insert into drug_regulons (drug_id, regulon_id) values (%s,%s)', [drug_id, regulon_id])
                seen.add(key)
        conn.commit()

def import_drug_programs(conn, df, drugs, programs):
    with conn.cursor() as cur:
        cur.execute('select count(*) from drug_programs')
        num_entries = cur.fetchone()[0]
        if num_entries > 0:
            return  # skip
        seen = set()
        for index, row in df.iterrows():
            drug_id = drugs[row['Drug']]
            program_id = programs['P-%d' % row['Programs']]
            key = '%d-%d' % (drug_id, program_id)
            if key not in seen:
                cur.execute('insert into drug_programs (drug_id,program_id) values (%s,%s)', [drug_id, program_id])
                seen.add(key)
        conn.commit()


def import_drug_data(conn, regulons, programs, genes, args):
    df = pd.read_csv(os.path.join(args.indir, 'Drugs_Mapped_to_Network_for_Portal.csv'),
                     sep=',', header=0)
    drugs = set(df['Drug'])
    drug_types = import_drug_types(conn, set(df['drugType']))
    action_types = import_action_types(conn, set(df['actionType']))
    mechanisms_of_action = import_mechanisms_of_action(conn, set(df["mechanismOfAction"]))
    target_type_models = import_target_type_models(conn, set(df["TargetTypeModel"]))

    drugs = import_drugs(conn, df, drug_types, action_types, mechanisms_of_action, target_type_models)
    import_drug_targets(conn, df, drugs, genes)
    import_drug_regulons(conn, df, drugs, regulons)
    import_drug_programs(conn, df, drugs, programs)


"""
def import_target_class_enrichment(conn, df, regulons):
    print('import target class enrichment')

    regulon2targetclass = defaultdict(set)
    all_targetclasses = set()
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_target_class')
        if cur.fetchone()[0] > 0:
            print('regulon target class enrichment found, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['TargetClassEnrichment']
            pval = row['TargetClass_p-value']
            try:
                if np.isnan(enrich):
                    enrich = None
            except:
                pass
            if enrich is not None:
                all_targetclasses.add(enrich)
                regulon2targetclass[regulon_id].add((enrich, pval))

        tc_ids = {}
        for tc in all_targetclasses:
            cur.execute('insert into target_classes (name) values (%s)', [tc])
            tc_ids[tc] = cur.lastrowid
        for regulon_id, target_classes in regulon2targetclass.items():
            for tc, pval in target_classes:
                tc_id = tc_ids[tc]
                cur.execute('insert into regulon_target_class (regulon_id,target_class_id,pval) values (%s,%s,%s)',
                            [regulon_id, tc_id, pval])
        conn.commit()
"""

"""
def import_hallmarks_enrichment(conn, df, regulons):
    print('import hallmarks enrichment')
    reg_hallmarks = defaultdict(set)
    reg_lin = defaultdict(set)
    all_hallmarks = set()
    all_lins = set()
    hallmark_ids = {}
    lin_hallmark_ids = {}

    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_hallmarks')
        if cur.fetchone()[0] > 0:
            print('entries exist, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['HallmarksEnrichment']
            lin_enrich = row['LinHallmarks']
            enrich = parse_list(enrich)
            lin_enrich = parse_list(lin_enrich)

            all_hallmarks.update(enrich)
            all_lins.update(lin_enrich)

            reg_hallmarks[regulon_id].update(enrich)
            reg_lin[regulon_id].update(lin_enrich)

        for hm in all_hallmarks:
            cur.execute('insert into hallmarks (name) values (%s)', [hm])
            hallmark_ids[hm] = cur.lastrowid

        for hm in all_lins:
            cur.execute('insert into lin_hallmarks (name) values (%s)', [hm])
            lin_hallmark_ids[hm] = cur.lastrowid

        for regulon_id, hallmarks in reg_hallmarks.items():
            for hm in hallmarks:
                hm_id = hallmark_ids[hm]
                cur.execute('insert into regulon_hallmarks (regulon_id,hallmark_id) values (%s,%s)',
                            [regulon_id, hm_id])

        for regulon_id, hallmarks in reg_lin.items():
            for hm in hallmarks:
                hm_id = lin_hallmark_ids[hm]
                cur.execute('insert into regulon_lin_hallmarks (regulon_id,hallmark_id) values (%s,%s)',
                            [regulon_id, hm_id])
        conn.commit()
"""

"""
def import_mechanism_of_action_enrichment(conn, df, regulons):
    print('import mechanism of action enrichment')

    all_moas = set()
    regulon_moas = defaultdict(set)
    moa_ids = {}

    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_mechanism_of_action')
        if cur.fetchone()[0] > 0:
            print('entries found, skipping')
            return

        for index, row in df.iterrows():
            regulon = row['Regulon']
            regulon_id = regulons[regulon][0]
            enrich = row['MechanismOfActionEnrichment']
            enrich = parse_list(enrich)
            all_moas.update(enrich)
            regulon_moas[regulon_id].update(enrich)

        for moa in all_moas:
            cur.execute('insert into mechanisms_of_action (name) values (%s)', [moa])
            moa_ids[moa] = cur.lastrowid

        for regulon_id, moas in regulon_moas.items():
            for moa in moas:
                moa_id = moa_ids[moa]
                cur.execute('insert into regulon_mechanism_of_action (regulon_id,mechanism_of_action_id) values (%s,%s)',
                            [regulon_id, moa_id])

        conn.commit()
"""

DBNAME = 'gbm_api'

HEADERS = [
    "Mutation", "Regulator", "Regulon", "MutationRegulatorEdge", "-log10(p)_MutationRegulatorEdge",
    "RegulatorRegulon_Spearman_R", "RegulatorRegulon_Spearman_p-value",
    "Regulon_stratification_t-statistic-log10(p)_Regulon_stratification	Fraction_of_edges_correctly_aligned",
    "TranscriptionalProgram", "HazardRatio", "HazardRatioPval", "RegulonGenes",
    "DrugEnrichment",
    "HallmarksEnrichment", "LinHallmarks",
    "TargetClassEnrichment", "TargetClass_p-value",
    "MechanismOfActionEnrichment", "MechanismOfAction_p-valuemirv_miRNA",
    "PITA_miRNA", "PITA_pval", "TargetScan_miRNA", "TargetScan_pval"
]

def dbconn():
    return MySQLdb.connect(host="localhost", user="admin", passwd="password",
                           db=DBNAME)


def read_synonyms(idconv_file):
    df = pd.read_csv(idconv_file, sep='\t')
    ens2pref = {}
    ens2entrez = {}
    for index, row in df.iterrows():
        ensembl = row['ensembl']
        try:
            entrez = int(row['entrez'].strip())
            ens2entrez[ensembl] = entrez
        except:
            pass
        try:
            preferred = row['preferred']
            ens2pref[ensembl] = preferred
        except:
            pass
    return ens2pref, ens2entrez



def read_cox_map(args):
    # make cox map
    cox_map = {}
    cox_map_df = pd.read_csv(os.path.join(args.indir, "CoxProportionalHazardsRegulons.csv"),
                             sep=',', header=0, index_col=0)
    for regulon, row in cox_map_df.iterrows():
        cox_map['R-%d' % regulon] = {'HazardRatio': row['HR'], 'p-value': row['p-value']}
    return cox_map

def read_regulons(args):
    # read the regulon map
    with open(os.path.join(args.indir, "regulons.json")) as infile:
        regulon_map = json.load(infile)
    # create a new regulon map, replacing the keys with "R-" prefixed
    new_regulon_map = {}
    for r, reg_genes in regulon_map.items():
        new_regulon_map['R-%s' % r] = reg_genes
    return new_regulon_map


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('indir', help="input directory")
    parser.add_argument('idconv', help="id conversion file")
    args = parser.parse_args()
    conn = dbconn()

    filtered_causals = os.path.join(args.indir, "filteredCausalResults.csv")
    df = pd.read_csv(filtered_causals, sep=',', header=0, index_col=0)

    cox_map = read_cox_map(args)
    ens2pref, ens2entrez = read_synonyms(args.idconv)
    regulon_map = read_regulons(args)

    gene_names = set()
    for regulon_genes in regulon_map.values():
        gene_names.update(regulon_genes)

    # Step 1: collect the Regulator field
    fill_roles(conn)
    tfs = import_tfs(conn, df)
    mutations = import_mutations(conn, df)
    genes = import_genes(conn, gene_names, ens2pref, ens2entrez)
    regulons = import_regulons(conn, regulon_map, cox_map, mutations)
    import_mutation_regulator(conn, df, regulons, mutations, tfs)
    import_regulon_regulator(conn, df, regulons, tfs)
    import_regulon_genes(conn, df, regulon_map, regulons, genes)
    programs = import_transcriptional_programs(conn, regulons, args)

    import_regulons_programs_genes_disease_mapping(conn,
                                                   regulons, programs, genes,
                                                   args)
    import_drug_data(conn, regulons, programs, genes, args)


    """
    import_drug_enrichment(conn, df, regulons)
    import_target_class_enrichment(conn, df, regulons)
    import_hallmarks_enrichment(conn, df, regulons)
    import_mechanism_of_action_enrichment(conn, df, regulons)
    """
