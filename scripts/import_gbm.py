#!/usr/bin/env python3

import os
import argparse
import MySQLdb
import pandas as pd
import numpy as np
import json
from collections import defaultdict


DESCRIPTION = "import_gbm - import filtered causal results to GBM backend"


def fill_roles(conn):
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_regulator_roles')
        # only insert if not exist
        if cur.fetchone()[0] == 0:
            cur.execute("insert into regulon_regulator_roles (id, name) values (1, 'activates')")
            cur.execute("insert into regulon_regulator_roles (id, name) values (2, 'represses')")
            cur.execute("insert into regulon_mutation_regulator_roles (id, name) values (1, 'down-regulates')")
            cur.execute("insert into regulon_mutation_regulator_roles (id, name) values (2, 'up-regulates')")
            conn.commit()

def import_regulators(conn, df, genes):
    """Import regulators"""
    regulators = set(df['Regulator'])
    result = {}
    with conn.cursor() as cur:
        cur.execute('select id,ensembl_id from genes where is_regulator=1')
        for pk, ensembl_id in cur.fetchall():
            result[ensembl_id] = pk
        for ens in regulators:
            if ens in genes:
                gene_id = genes[ens][0]
                cur.execute('update genes set is_regulator=1 where id=%s', [gene_id])
                result[ens] = gene_id
            else:
                cur.execute("insert into genes (ensembl_id,is_regulator) values (%s,1)", [ens])
                regulator_id = cur.lastrowid
                genes[ens] = (regulator_id, None, None)
                result[ens] = regulator_id
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


def import_genes(conn, ens_genes, ens2pref, ens2entrez):
    """
    genes: id: int, ensembl_id: string, entrez_id: string, preferred: string
    In the import document this is the Regulator
    """
    print("# unique genes found: %d" % len(ens_genes))

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
            for ens in sorted(ens_genes):
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
    """regulons: id, name, cox_hazard_ratio"""
    result = {}
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulons')
        if cur.fetchone()[0] > 0:
            print('Regulons found, reading from database')
            cur.execute('select id,name,cox_hazard_ratio from regulons')
            for pk, regulon, hr in cur.fetchall():
                result[regulon] = (pk, hr)
        else:
            print('Regulons not found, inserting into database')
            for regulon in regulon_map.keys():
                cox_hr = cox_map[regulon]['HazardRatio']
                cur.execute('select count(*) from regulons where name=%s', [regulon])
                if cur.fetchone()[0] == 0:
                    cur.execute('insert into regulons (name,cox_hazard_ratio) values (%s,%s)',
                                [regulon, cox_hr])
                    result[regulon] = (cur.lastrowid, cox_hr)
            conn.commit()
    return result


def import_mutation_regulator(conn, df, regulons, mutations, tfs):
    """table: regulon_mutation_regulator, field: MutationRegulatorEdge
    id, regulon_id, mutation_id, regulator_id, role_id"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_mutation_regulator')
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
                    role_id = 1  # down-regulates
                else:
                    role_id = 2  # up-regulates
                cur.execute('insert into regulon_mutation_regulator (regulon_id,mutation_id,regulator_id,role_id) values (%s,%s,%s,%s)',
                            [regulon_id, mutation_id, regulator_id, role_id])
            conn.commit()
        else:
            print('skip inserting mutation regulator edges')


def import_regulon_regulator(conn, df, regulons, tfs):
    """table: bc_tf (regulon_id,tf_id,role)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_regulator')
        if cur.fetchone()[0] == 0:
            print('insert regulon-regulator into database')
            for index, row in df.iterrows():
                regulon = row['Regulon']
                regulator = row['Regulator']
                regulon_id = regulons[regulon][0]
                regulator_id = tfs[regulator]
                spearman_r = row['RegulatorRegulon_Spearman_R']
                if spearman_r > 0:
                    role_id = 1  # 1 activates
                else:
                    role_id = 2  # 2 represses
                #print(row['RegulatorRegulon_Spearman_p-value'])
                cur.execute('select count(*) from regulon_regulator where regulon_id=%s and regulator_id=%s',
                            [regulon_id, regulator_id])
                if cur.fetchone()[0] == 0:
                    # only insert once
                    cur.execute('insert into regulon_regulator (regulon_id,regulator_id,role_id) values (%s,%s,%s)',
                                [regulon_id, regulator_id, role_id])
            conn.commit()
        else:
            print('regulon regulator relations exist, skip')


def import_regulon_genes(conn, df, regulon_map, regulons, genes):
    """table: regulon_genes (regulon_id, gene_id)"""
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_genes')
        if cur.fetchone()[0] == 0:
            print('insert regulon genes into database')
            for regulon, regulon_genes in regulon_map.items():
                regulon_id = regulons[regulon][0]
                for gene in regulon_genes:
                    gene_id = genes[gene][0]
                    cur.execute('select count(*) from regulon_genes where regulon_id=%s and gene_id=%s',
                                [regulon_id, gene_id])
                    if cur.fetchone()[0] == 0:  # prevent double insertions in the same transaction
                        cur.execute('insert into regulon_genes(regulon_id,gene_id) values (%s,%s)',
                                    [regulon_id, gene_id])
            conn.commit()
        else:
            print('regulon genes exist, skip')


def import_transcriptional_programs(conn, regulons, args):
    with open(os.path.join(args.indir, 'transcriptional_programs.json')) as infile:
        program2regulons = json.load(infile)
    print('import transcriptional programs')
    result = {}
    programs = {}
    with conn.cursor() as cur:
        for program_number, prog_regulons in program2regulons.items():
            # add to database if necessary
            progname = "P-%s" % program_number
            if progname not in programs:
                cur.execute('insert into trans_programs (name) values (%s)', [progname])
                programs[progname] = cur.lastrowid
            for regulon_number in prog_regulons:
                regulon = "R-%s" % regulon_number
                regulon_id = regulons[regulon][0]
                cur.execute('insert into regulon_programs (regulon_id, program_id) values (%s,%s)', [regulon_id, programs[progname]])
        conn.commit()
    return programs

def import_regulons_programs_genes_disease_mapping(conn,
                                                   regulons, programs, genes,
                                                   args):
    df = pd.read_csv(os.path.join(args.indir, 'regulons_programs_genes_disease_mapping.csv'),
                     sep=',', header=0)
    with conn.cursor() as cur:
        cur.execute('select count(*) from regulon_program_genes')
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
            cur.execute('insert into regulon_program_genes (regulon_id,program_id,gene_id,is_disease_relevant) values (%s,%s,%s,%s)', [regulon_id, program_id, gene_id, 1 if is_disease_relevant != 'False' else 0])
        conn.commit()


def read_catalog_table_as_map(conn, table):
    result = {}
    with conn.cursor() as cur:
        cur.execute('select id, name from %s' % table)
        for pk, name in cur.fetchall():
            result[name] = pk
        return result

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


def import_cmflows(conn, regulons, mutations, genes, tfs, filename, has_pathway=False):
    with conn.cursor() as cur:
        cmflow_types = read_catalog_table_as_map(conn, "cm_flow_types")
        cmf_pathways = read_catalog_table_as_map(conn, "cmf_pathways")
        symbol2gene_id = {}
        cur.execute('select id,preferred from genes')
        for pk, symbol in cur.fetchall():
            symbol2gene_id[symbol] = pk

        df = pd.read_csv(os.path.join(args.indir, filename),
                         sep=',', header=0)
        # 1. import cm_flow_types, note: load from database first
        df_flow_types = set(df['CMFlowType'])
        for ftype in df_flow_types:
            if not ftype in cmflow_types:
                cur.execute('insert into cm_flow_types (name) values (%s)', [ftype])
                pk = cur.lastrowid
                cmflow_types[ftype] = pk
        conn.commit()

        # 2. import cmf_pathways, note: load from database
        if has_pathway:
            df_pathways = set(df['Pathway'])
            for pway in df_pathways:
                if not pway in cmf_pathways:
                    cur.execute('insert into cmf_pathways (name) values (%s)', [pway])
                    pk = cur.lastrowid
                    cmf_pathways[pway] = pk
            conn.commit()

        # 3. import mutations, note: make sure they exist
        df_mutations = set(df['Mutation'])
        for mutation in df_mutations:
            if mutation not in mutations:
                cur.execute('insert into mutations (name) values (%s)', [mutation])
                pk = cur.lastrowid
                mutations[mutation] = pk

        for index, row in df.iterrows():
            cmf_id = row['CMF_ID']
            cmf_name = row['CMFlow']
            cmf_type_id = cmflow_types[row['CMFlowType']]
            mutation_id = mutations[row['Mutation']]
            if has_pathway:
                cmf_pathway_id = cmf_pathways[row['Pathway']]
            else:
                cmf_pathway_id = None

            # import mutation genes
            # Brutal truth: Sometimes, only the symbol exists and not the
            # EnsEMBL ID, so we need to add that
            mut_ens = row['MutationEnsembl']
            mut_sym = row['MutationSymbol']
            if not has_pathway:
                if type(mut_ens) != str and np.isnan(mut_ens):
                    #print("mutation gene is NaN for symbol: '%s'" % mut_sym)
                    if mut_sym not in symbol2gene_id:
                        cur.execute('insert into genes (preferred,is_mutation) values (%s,1)', [mut_sym])
                        mutation_gene_id = cur.lastrowid
                        genes[mut_ens] = (mutation_gene_id, None, mut_sym)
                        symbol2gene_id[mut_sym] = mutation_gene_id
                    else:
                        mutation_gene_id = symbol2gene_id[mut_sym]
                else:
                    if mut_ens not in genes:
                        # Insert as a new gene
                        cur.execute('insert into genes (ensembl_id,preferred,is_mutation) values (%s,%s,1)',
                                    [mut_ens, mut_sym])
                        mutation_gene_id = cur.lastrowid
                        genes[mut_ens] = (mutation_gene_id, None, mut_sym)
                    else:
                        mutation_gene_id = genes[mut_ens][0]
                        cur.execute('update genes set is_mutation=1 where id=%s', [mutation_gene_id])
                    conn.commit()
            else:
                # Pathway file does not have mutation gene info
                mutation_gene_id = None

            regulator_id = genes[row['RegulatorEnsembl']][0]
            regulon_id = regulons["R-%d" % row['Regulon']][0]
            bc_mutation_tf_role_id = 2 if row['MutationRegulatorEdge'] == 1 else 1
            bc_mutation_tf_pvalue = row['-log10(p)_MutationRegulatorEdge']
            bc_tf_spearman_r = row['RegulatorRegulon_Spearman_R']
            bc_tf_spearman_pvalue = row['RegulatorRegulon_Spearman_p-value']
            bc_t_statistic = row['Regulon_stratification_t-statistic']
            bc_log10_p_stratification = row['-log10(p)_Regulon_stratification']
            fraction_edges_correctly_aligned = row['Fraction_of_edges_correctly_aligned']
            fraction_aligned_diffexp_edges = row['Fraction_of_aligned_and_diff_exp_edges']
            num_downstream_regulons = row['number_downstream_regulons']
            num_diffexp_regulons = row['number_differentially_expressed_regulons']

            query = """insert into cm_flows (cmf_id,cmf_name,cmf_type_id,
            mutation_id,cmf_pathway_id,mutation_gene_id,regulator_id,regulon_id,
            regulon_mutation_regulator_role_id,regulon_mutation_regulator_pvalue,
            regulon_regulator_spearman_r,regulon_regulator_spearman_pvalue,
            regulon_t_statistic,regulon_log10_p_stratification,
            fraction_edges_correctly_aligned,fraction_aligned_diffexp_edges,
            num_downstream_regulons,num_diffexp_regulons) values
            (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"""
            cur.execute(query, [cmf_id, cmf_name, cmf_type_id, mutation_id, cmf_pathway_id,
                                mutation_gene_id, regulator_id, regulon_id,
                                bc_mutation_tf_role_id, bc_mutation_tf_pvalue,
                                bc_tf_spearman_r, bc_tf_spearman_pvalue,
                                bc_t_statistic, bc_log10_p_stratification,
                                fraction_edges_correctly_aligned, fraction_aligned_diffexp_edges,
                                num_downstream_regulons, num_diffexp_regulons])


DBNAME = 'gbm_api'


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
    mutations = import_mutations(conn, df)
    genes = import_genes(conn, gene_names, ens2pref, ens2entrez)

    regulators = import_regulators(conn, df, genes)
    regulons = import_regulons(conn, regulon_map, cox_map, mutations)
    import_mutation_regulator(conn, df, regulons, mutations, regulators)
    import_regulon_regulator(conn, df, regulons, regulators)
    import_regulon_genes(conn, df, regulon_map, regulons, genes)
    programs = import_transcriptional_programs(conn, regulons, args)

    import_regulons_programs_genes_disease_mapping(conn,
                                                   regulons, programs, genes,
                                                   args)
    import_drug_data(conn, regulons, programs, genes, args)
    import_cmflows(conn, regulons, mutations, genes, regulators, "diseaseRelevantCMFlowsGenes.csv", False)
    import_cmflows(conn, regulons, mutations, genes, regulators, "diseaseRelevantCMFlowsPathways.csv", True)
